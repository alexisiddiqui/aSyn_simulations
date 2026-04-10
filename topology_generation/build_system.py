"""Build GROMACS simulation boxes for each condition.

Pipeline per condition:
  1. assemble_topology   — combine protein.top + ligand ITPs
  2. gmx editconf        — set box size/shape, center protein
  3. gmx insert-molecules — insert buffer ligands
  4. gmx solvate         — fill with OPC water
  5. gmx grompp+genion   — add concentration ions (one pass per species)
  6. gmx grompp+genion   — neutralize with Cl-/Na+
"""

import argparse
import os
from pathlib import Path

import yaml

from system_builder_utils import (
    assemble_topology,
    compute_box_volume,
    compute_buffer_counts,
    compute_ion_counts,
    run_editconf,
    run_genion,
    run_grompp_for_genion,
    run_insert_molecules,
    run_solvate,
)
from topology_utils import load_config


def process_condition(
    condition_name: str,
    topo_cfg: dict,
    sys_cfg: dict,
    output_root: Path,
    topo_output_root: Path,
    repo_root: Path,
    gmx_path: str,
    ions_mdp: Path,
    force: bool = False,
) -> Path:
    """Build a solvated, ionised simulation box for one condition.

    Args:
        condition_name: Condition key (e.g. 'tris', 'extracellular')
        topo_cfg: topology_config.yaml parsed dict
        sys_cfg: system_config.yaml parsed dict
        output_root: system_output/ directory
        topo_output_root: topology_output/ directory
        repo_root: Repository root
        gmx_path: Path to gmx executable
        ions_mdp: Path to config/GROMACS/ions.mdp
        force: If True, re-run even if outputs already exist

    Returns:
        Path to condition output directory
    """
    work_dir = output_root / condition_name
    work_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Building system: {condition_name}")
    print(f"{'='*60}")

    # Check if already complete
    final_gro = work_dir / "system.gro"
    final_top = work_dir / "system.top"
    if final_gro.exists() and final_top.exists() and not force:
        print(f"  Output already exists. Skipping (use --force to re-run)")
        return work_dir

    # --- Resolve inputs from topology_config ---
    cond_topo = topo_cfg["conditions"][condition_name]
    ph = cond_topo["ph"]

    protein_dir = topo_output_root / "protein" / f"ph_{ph}"
    protein_gro = protein_dir / "protein.gro"
    protein_top = protein_dir / "protein.top"
    protein_posre = protein_dir / "posre.itp"

    for p in [protein_gro, protein_top]:
        if not p.exists():
            raise FileNotFoundError(
                f"Protein topology not found for pH {ph}: {p}\n"
                f"Run generate_protein_topology.py first."
            )

    # --- Resolve ligand inputs ---
    ligands_in_condition = cond_topo.get("ligands", [])
    ligand_info = []  # (resname, itp_path, top_path)
    for resname in ligands_in_condition:
        ligand_dir = topo_output_root / "ligands" / resname
        itp_path = ligand_dir / f"{resname}_GMX.itp"
        top_path = ligand_dir / f"{resname}_GMX.top"
        for p in [itp_path, top_path]:
            if not p.exists():
                raise FileNotFoundError(
                    f"Ligand topology not found: {p}\n"
                    f"Run generate_ligand_topology.py first."
                )
        ligand_info.append((resname, itp_path, top_path))

    # --- Box geometry and counts ---
    box_size = sys_cfg["box"]["size_nm"]
    box_shape = sys_cfg["box"]["shape"]
    volume = compute_box_volume(box_size, box_shape)

    buffer_conc = sys_cfg["buffer"]["concentration_mm"]
    sys_cond = sys_cfg["conditions"][condition_name]
    buffer_species = sys_cond.get("buffer_species", {})
    ions_cfg = sys_cond.get("ions", {})

    buffer_counts = compute_buffer_counts(buffer_species, buffer_conc, volume)
    ion_list = compute_ion_counts(ions_cfg, volume)

    neut_ion = sys_cfg["neutralization"]["ion"]
    neut_charge = sys_cfg["neutralization"]["charge"]

    print(f"\n  pH: {ph}")
    print(f"  Box: {box_size} nm {box_shape} (V = {volume:.1f} nm³)")
    if buffer_counts:
        print(f"  Buffer molecules: {buffer_counts}")
    else:
        print(f"  Buffer molecules: none")
    if ion_list:
        for ion, count, charge in ion_list:
            print(f"  Ion {ion} (charge {charge:+d}): {count} molecules")
    else:
        print(f"  Ions: none (neutralization only)")

    # --- Step 1: Assemble topology ---
    print(f"\n[1] Assembling topology...")
    system_top = assemble_topology(
        protein_top=protein_top,
        ligand_info=ligand_info,
        buffer_counts=buffer_counts,
        work_dir=work_dir,
    )
    print(f"  -> {system_top}")

    # --- Step 2: editconf ---
    print(f"\n[2] Defining simulation box (editconf)...")
    box_gro = work_dir / "box.gro"
    run_editconf(gmx_path, protein_gro, box_gro, box_size, box_shape, work_dir)
    print(f"  -> {box_gro}")

    # --- Step 3: Insert buffer ligands ---
    current_gro = box_gro
    if buffer_counts:
        print(f"\n[3] Inserting buffer ligands...")
        for i, (resname, itp_path, _) in enumerate(ligand_info):
            count = buffer_counts.get(resname, 0)
            if count == 0:
                continue

            ligand_gro = (
                topo_output_root / "ligands" / resname
                / f"{resname}.amb2gmx" / f"{resname}_GMX.gro"
            )
            if not ligand_gro.exists():
                raise FileNotFoundError(f"Ligand .gro not found: {ligand_gro}")

            out_gro = work_dir / f"after_{resname}.gro"
            label = f"insert_{resname}"
            print(f"  Inserting {count} x {resname}...")
            actual = run_insert_molecules(
                gmx_path, current_gro, ligand_gro, count, out_gro, work_dir, label
            )
            print(f"  -> {out_gro} ({actual} inserted)")
            current_gro = out_gro
    else:
        print(f"\n[3] No buffer ligands to insert.")

    # --- Step 4: Solvate ---
    print(f"\n[4] Solvating system...")
    solvated_gro = work_dir / "solvated.gro"
    water_model = topo_cfg["protein"]["water_model"]
    run_solvate(gmx_path, current_gro, system_top, solvated_gro, work_dir, water_model)
    print(f"  -> {solvated_gro}")
    current_gro = solvated_gro

    # --- Step 5: Add concentration ions (one genion pass per species) ---
    if ion_list:
        print(f"\n[5] Adding concentration ions...")
        for ion_name, n, charge in ion_list:
            label = f"genion_{ion_name}"
            tpr = work_dir / "ions.tpr"
            out_gro = work_dir / f"after_{ion_name}.gro"

            print(f"  Adding {n} x {ion_name} (charge {charge:+d})...")

            if charge > 0:
                pname, pcharge = ion_name, charge
                nname, ncharge_val = neut_ion, neut_charge
                n_pos, n_neg = n, 0
            else:
                pname, pcharge = "NA", 1
                nname, ncharge_val = ion_name, charge
                n_pos, n_neg = 0, n

            run_grompp_for_genion(
                gmx_path, ions_mdp, current_gro, system_top, tpr, work_dir,
                f"grompp_{label}"
            )
            run_genion(
                gmx_path, tpr, system_top, out_gro, work_dir, label,
                n_pos=n_pos, pname=pname, pcharge=pcharge,
                n_neg=n_neg, nname=nname, ncharge=ncharge_val,
            )
            print(f"  -> {out_gro}")
            current_gro = out_gro
    else:
        print(f"\n[5] No concentration ions to add.")

    # --- Step 6: Neutralization ---
    print(f"\n[6] Neutralizing system...")
    tpr = work_dir / "ions_neutral.tpr"
    run_grompp_for_genion(
        gmx_path, ions_mdp, current_gro, system_top, tpr, work_dir,
        "grompp_neutral"
    )
    run_genion(
        gmx_path, tpr, system_top, final_gro, work_dir, "genion_neutral",
        pname="NA", pcharge=1,
        nname=neut_ion, ncharge=neut_charge,
        neutral=True,
    )
    print(f"  -> {final_gro}")

    # --- Validation ---
    print(f"\n  Validating outputs...")
    validate_system_outputs(work_dir, condition_name)

    print(f"\n  Condition {condition_name!r} complete")
    return work_dir


def validate_system_outputs(work_dir: Path, condition_name: str) -> None:
    """Basic validation of system build outputs.

    Args:
        work_dir: Condition output directory
        condition_name: For error messages
    """
    gro = work_dir / "system.gro"
    top = work_dir / "system.top"

    if not gro.exists():
        raise FileNotFoundError(f"system.gro not found in {work_dir}")
    if not top.exists():
        raise FileNotFoundError(f"system.top not found in {work_dir}")

    # Check atom count in .gro
    lines = gro.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"system.gro appears empty")
    try:
        n_atoms = int(lines[1].strip())
    except ValueError:
        raise ValueError(f"Cannot parse atom count from system.gro")

    print(f"  system.gro: {n_atoms} atoms")

    # Check for errors in any log file
    for log in work_dir.glob("*.log"):
        content = log.read_text()
        for line in content.splitlines():
            if "Fatal error" in line or "fatal error" in line:
                print(f"  WARNING: potential error in {log.name}: {line.strip()}")

    # Check neutralization log for charge info
    neut_log = work_dir / "genion_neutral.log"
    if neut_log.exists():
        content = neut_log.read_text()
        for line in content.splitlines():
            if "charge" in line.lower() and ("system" in line.lower() or "total" in line.lower()):
                print(f"  {line.strip()}")

    print(f"  Validation passed")


def main():
    parser = argparse.ArgumentParser(
        description="Build GROMACS simulation boxes for each condition"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/topology_config.yaml"),
        help="Path to topology_config.yaml",
    )
    parser.add_argument(
        "--sys-config",
        type=Path,
        default=Path("config/system_config.yaml"),
        help="Path to system_config.yaml",
    )
    parser.add_argument(
        "--condition",
        type=str,
        default=None,
        help="Build a single condition (default: all conditions)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="System output folder (default: system_output/ relative to repo root)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if outputs already exist",
    )
    args = parser.parse_args()

    # Load configs
    topo_cfg = load_config(args.config)
    repo_root = Path(topo_cfg["paths"]["repo_root"])
    os.chdir(repo_root)

    with open(args.sys_config) as f:
        sys_cfg = yaml.safe_load(f)

    topo_output_root = repo_root / topo_cfg["paths"]["output_dir"]
    output_root = args.output if args.output else repo_root / "system_output"
    output_root.mkdir(parents=True, exist_ok=True)

    gmx_path = topo_cfg["paths"]["gmx"]
    if not Path(gmx_path).exists():
        raise FileNotFoundError(f"GROMACS executable not found: {gmx_path}")

    ions_mdp = repo_root / "config" / "GROMACS" / "ions.mdp"
    if not ions_mdp.exists():
        raise FileNotFoundError(f"ions.mdp not found: {ions_mdp}")

    # Determine which conditions to build
    all_conditions = list(sys_cfg["conditions"].keys())
    if args.condition:
        if args.condition not in all_conditions:
            raise ValueError(
                f"Unknown condition: {args.condition!r}. "
                f"Available: {all_conditions}"
            )
        conditions = [args.condition]
    else:
        conditions = all_conditions

    box_size = sys_cfg["box"]["size_nm"]
    box_shape = sys_cfg["box"]["shape"]

    print(f"\nSystem Building Pipeline")
    print(f"Repository:  {repo_root}")
    print(f"Topology in: {topo_output_root}")
    print(f"Output:      {output_root}")
    print(f"Box:         {box_size} nm {box_shape}")
    print(f"Conditions:  {conditions}")

    for condition_name in conditions:
        process_condition(
            condition_name=condition_name,
            topo_cfg=topo_cfg,
            sys_cfg=sys_cfg,
            output_root=output_root,
            topo_output_root=topo_output_root,
            repo_root=repo_root,
            gmx_path=gmx_path,
            ions_mdp=ions_mdp,
            force=args.force,
        )

    print(f"\n{'='*60}")
    print(f"All conditions complete")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
