"""Run GROMACS simulation pipeline for aSyn MD.

Pipeline stages: energy minimisation → equilibration → relaxation → production
Per condition: 5 conditions × (1 minim + 5 equil + 1 relax × 3 replicates
              + 1 smoke test × 3 + 1 production × 3) = 5 × 20 = 100 stages total

Output: simulation_output/{condition}/{stage}/md.{tpr,gro,cpt,xtc,edr,log}
"""

import argparse
import os
from pathlib import Path
from dataclasses import asdict

import yaml

from simulation_runner_utils import (
    StageConfig,
    run_make_ndx,
    run_grompp,
    run_mdrun,
)
from topology_utils import load_config


def load_stages(stages_config_path: Path) -> list[StageConfig]:
    """Load stage definitions from YAML.

    Returns:
        List of StageConfig objects in order.
    """
    with open(stages_config_path) as f:
        data = yaml.safe_load(f)

    stages = []
    for stage_dict in data['stages']:
        stages.append(StageConfig(**stage_dict))

    return stages


def resolve_posres_itp(repo_root: Path, condition_name: str, topo_cfg: dict, work_dir: Path) -> Path:
    """Copy posre.itp from topology_output to work_dir.

    Args:
        repo_root: Repository root
        condition_name: Condition name (e.g. 'tris')
        topo_cfg: Parsed topology_config.yaml
        work_dir: Simulation work directory

    Returns:
        Path to posre.itp in work_dir
    """
    cond_topo = topo_cfg['conditions'][condition_name]
    ph = cond_topo['ph']

    topo_output_root = repo_root / topo_cfg['paths']['output_dir']
    posre_src = topo_output_root / 'protein' / f'ph_{ph}' / 'posre.itp'

    if not posre_src.exists():
        raise FileNotFoundError(f"posre.itp not found: {posre_src}")

    posre_dst = work_dir / 'posre.itp'
    if not posre_dst.exists():
        posre_dst.write_text(posre_src.read_text())
        print(f"  Copied posre.itp → {posre_dst.name}")

    return posre_dst


def process_condition(
    condition_name: str,
    topo_cfg: dict,
    stages: list[StageConfig],
    repo_root: Path,
    gmx_path: str,
    system_root: Path,
    output_root: Path,
    n_replicates: int = 3,
    through_stage: str = None,
    single_stage: str = None,
    force: bool = False,
) -> None:
    """Build and run all simulation stages for one condition.

    Args:
        condition_name: Condition key (e.g. 'tris')
        topo_cfg: Parsed topology_config.yaml
        stages: List of StageConfig objects
        repo_root: Repository root
        gmx_path: Path to gmx executable
        stages_config: Path to simulation_stages.yaml (for resolving MDP paths)
        output_root: simulation_output/ directory
        n_replicates: Number of replicates for replicated stages
        through_stage: Stop after this stage (e.g. '6_equil')
        single_stage: Run only this stage (implies force=True for it)
        force: Re-run all stages even if outputs exist
    """
    work_dir = output_root / condition_name
    work_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Condition: {condition_name}")
    print(f"{'='*60}")

    # Copy system files
    system_src = system_root / condition_name
    system_gro = system_src / 'system.gro'
    system_top = system_src / 'system.top'

    if not system_gro.exists() or not system_top.exists():
        raise FileNotFoundError(
            f"System files not found for {condition_name}. "
            f"Run build_system.py first."
        )

    work_system_gro = work_dir / 'system.gro'
    work_system_top = work_dir / 'system.top'

    if not work_system_gro.exists():
        work_system_gro.write_text(system_gro.read_text())
        print(f"  Copied system.gro")

    if not work_system_top.exists():
        work_system_top.write_text(system_top.read_text())
        print(f"  Copied system.top")

    # Copy and resolve posre.itp
    resolve_posres_itp(repo_root, condition_name, topo_cfg, work_dir)

    # Generate index file
    ndx_file = work_dir / 'system.ndx'
    print(f"\n[Setup] Creating index file...")
    run_make_ndx(gmx_path, work_system_gro, ndx_file, work_dir)

    # Run stages
    print(f"\n[Stages]")

    stages_to_run = stages

    if single_stage:
        stages_to_run = [s for s in stages if s.name == single_stage]
        if not stages_to_run:
            raise ValueError(f"Unknown stage: {single_stage}")
        force = True
    elif through_stage:
        stage_names = [s.name for s in stages]
        if through_stage not in stage_names:
            raise ValueError(
                f"Unknown stage: {through_stage}. "
                f"Available: {stage_names}"
            )
        through_idx = stage_names.index(through_stage)
        stages_to_run = stages[:through_idx + 1]

    for stage in stages_to_run:
        if stage.replicated:
            run_replicated_stage(
                stage,
                work_dir,
                gmx_path,
                stages,
                n_replicates,
                force
            )
        else:
            run_single_stage(
                stage,
                work_dir,
                gmx_path,
                stages,
                force
            )

    print(f"\n{'='*60}")
    print(f"Condition {condition_name} complete")
    print(f"{'='*60}")


def run_single_stage(
    stage: StageConfig,
    work_dir: Path,
    gmx_path: str,
    stages_list: list[StageConfig],
    force: bool = False,
) -> Path:
    """Run a non-replicated stage.

    Returns:
        Path to stage output directory
    """
    stage_dir = work_dir / stage.name
    stage_dir.mkdir(parents=True, exist_ok=True)

    md_gro = stage_dir / 'md.gro'
    if md_gro.exists() and not force:
        print(f"  [{stage.name}] Already complete, skipping")
        return stage_dir

    print(f"  [{stage.name}]")

    # Resolve previous stage outputs
    if stage.name == '1_minim':
        prev_gro = work_dir / 'system.gro'
        prev_cpt = None
    else:
        # Find previous non-replicated stage in order
        stage_idx = next(i for i, s in enumerate(stages_list) if s.name == stage.name)
        prev_non_repl = [s for s in stages_list[:stage_idx] if not s.replicated]
        if not prev_non_repl:
            raise RuntimeError(f"Cannot determine previous stage for {stage.name}")
        prev_stage_name = prev_non_repl[-1].name
        prev_gro = work_dir / prev_stage_name / 'md.gro'
        prev_cpt = work_dir / prev_stage_name / 'md.cpt'

    if not prev_gro.exists():
        raise FileNotFoundError(f"Previous .gro not found: {prev_gro}")

    # Resolve MDP path
    mdp_path = Path(stage.mdp)
    if not mdp_path.is_absolute():
        mdp_path = Path(os.getcwd()) / mdp_path

    # grompp
    run_grompp(
        gmx_path,
        mdp_path,
        prev_gro,
        work_dir / 'system.top',
        work_dir / 'system.ndx',
        stage_dir / 'md.tpr',
        stage_dir,
        f"grompp_{stage.name}",
        ref_file=prev_gro if stage.use_posres_ref else None,
        cpt_file=prev_cpt if stage.continuation else None,
    )

    # mdrun
    run_mdrun(gmx_path, stage_dir / 'md.tpr', stage_dir, stage.name)

    return stage_dir


def run_replicated_stage(
    stage: StageConfig,
    work_dir: Path,
    gmx_path: str,
    stages_config: list[StageConfig],
    n_replicates: int,
    force: bool = False,
) -> None:
    """Run a replicated stage n_replicates times.

    Args:
        stage: StageConfig for this stage
        work_dir: Condition work directory
        gmx_path: Path to gmx executable
        stages_config: Full list of stages (to look up dependencies)
        n_replicates: Number of replicates
        force: Re-run even if outputs exist
    """
    print(f"  [{stage.name}] (×{n_replicates} replicates)")

    for rep in range(1, n_replicates + 1):
        stage_dir = work_dir / f"{stage.name}_rep{rep}"
        stage_dir.mkdir(parents=True, exist_ok=True)

        md_gro = stage_dir / 'md.gro'
        if md_gro.exists() and not force:
            print(f"    rep{rep}: Already complete, skipping")
            continue

        # Determine input: depends_on or previous non-replicated stage
        if stage.depends_on:
            prev_stage_name = stage.depends_on
            prev_dir = work_dir / f"{prev_stage_name}_rep{rep}"
        else:
            # Find last non-replicated stage before this one
            preceding = [s for s in stages_config if s.name < stage.name and not s.replicated]
            if not preceding:
                raise RuntimeError(f"No non-replicated stage before {stage.name}")
            prev_stage_name = preceding[-1].name
            prev_dir = work_dir / prev_stage_name

        prev_gro = prev_dir / 'md.gro'
        prev_cpt = prev_dir / 'md.cpt'

        if not prev_gro.exists():
            raise FileNotFoundError(f"Previous .gro not found: {prev_gro}")

        # Resolve MDP path
        mdp_path = Path(stage.mdp)
        if not mdp_path.is_absolute():
            mdp_path = Path(os.getcwd()) / mdp_path

        # grompp
        run_grompp(
            gmx_path,
            mdp_path,
            prev_gro,
            work_dir / 'system.top',
            work_dir / 'system.ndx',
            stage_dir / 'md.tpr',
            stage_dir,
            f"grompp_{stage.name}_rep{rep}",
            ref_file=prev_gro if stage.use_posres_ref else None,
            cpt_file=prev_cpt if stage.continuation else None,
        )

        # mdrun
        run_mdrun(gmx_path, stage_dir / 'md.tpr', stage_dir, f"{stage.name}_rep{rep}")

        print(f"    rep{rep}: complete")


def main():
    parser = argparse.ArgumentParser(
        description="Run GROMACS simulation pipeline for aSyn MD"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/topology_config.yaml"),
        help="Path to topology_config.yaml",
    )
    parser.add_argument(
        "--stages-config",
        type=Path,
        default=Path("config/simulation_stages.yaml"),
        help="Path to simulation_stages.yaml",
    )
    parser.add_argument(
        "--condition",
        type=str,
        default=None,
        help="Run single condition (default: all)",
    )
    parser.add_argument(
        "--through",
        type=str,
        default=None,
        help="Run up to and including this stage",
    )
    parser.add_argument(
        "--stage",
        type=str,
        default=None,
        help="Run only this stage (implies --force)",
    )
    parser.add_argument(
        "--n-replicates",
        type=int,
        default=3,
        help="Number of replicates for replicated stages",
    )
    parser.add_argument(
        "--system-output",
        type=Path,
        default=None,
        help="Folder containing built systems (default: system_output/ relative to repo root)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Simulation output folder (default: simulations/ inside --system-output)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if outputs already exist",
    )
    args = parser.parse_args()

    # Load configs
    topo_cfg = load_config(args.config)
    repo_root = Path(topo_cfg['paths']['repo_root'])
    os.chdir(repo_root)

    stages = load_stages(args.stages_config)

    gmx_path = topo_cfg['paths']['gmx']
    if not Path(gmx_path).exists():
        raise FileNotFoundError(f"GROMACS executable not found: {gmx_path}")

    system_root = args.system_output if args.system_output else repo_root / "system_output"
    output_root = args.output if args.output else system_root / "simulations"
    output_root.mkdir(parents=True, exist_ok=True)

    # Determine which conditions to run
    all_conditions = list(topo_cfg['conditions'].keys())
    if args.condition:
        if args.condition not in all_conditions:
            raise ValueError(
                f"Unknown condition: {args.condition}. "
                f"Available: {all_conditions}"
            )
        conditions = [args.condition]
    else:
        conditions = all_conditions

    print(f"\nSimulation Pipeline")
    print(f"Repository:  {repo_root}")
    print(f"Systems in:  {system_root}")
    print(f"Output:      {output_root}")
    print(f"Conditions:  {conditions}")
    print(f"Replicates:  {args.n_replicates}")

    for condition_name in conditions:
        process_condition(
            condition_name=condition_name,
            topo_cfg=topo_cfg,
            stages=stages,
            repo_root=repo_root,
            gmx_path=gmx_path,
            system_root=system_root,
            output_root=output_root,
            n_replicates=args.n_replicates,
            through_stage=args.through,
            single_stage=args.stage,
            force=args.force,
        )

    print(f"\n{'='*60}")
    print(f"All conditions complete")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
