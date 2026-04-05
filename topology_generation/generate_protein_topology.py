"""Generate GROMACS topologies for protein using pdb2gmx at different pH values."""

import os
import subprocess
from pathlib import Path

from topology_utils import load_config, preprocess_pdb_for_ph, run_cmd, validate_files_exist


def process_ph(
    ph: float,
    protein_cfg: dict,
    ph_prot_cfg: dict,
    output_root: Path,
    gmx_path: str,
    force: bool = False,
) -> Path:
    """Run pdb2gmx for one pH value.

    Args:
        ph: pH value (e.g., 7.4, 7.2, 4.9)
        protein_cfg: Config dict with 'input_pdb', 'forcefield', 'water_model'
        ph_prot_cfg: pH-specific protonation config ('his_method', 'overrides')
        output_root: Root output directory (topology_output)
        gmx_path: Path to gmx executable
        force: If True, re-run even if outputs exist

    Returns:
        Path to protein output directory for this pH
    """
    # pH as string for directory naming (handle floats)
    ph_str = str(ph).replace(".", "_")
    work_dir = output_root / "protein" / f"ph_{ph}"
    work_dir.mkdir(parents=True, exist_ok=True)

    pdb_src = protein_cfg["input_pdb"]
    forcefield = protein_cfg["forcefield"]
    water_model = protein_cfg["water_model"]

    print(f"\n{'='*60}")
    print(f"Processing protein at pH {ph}")
    print(f"Forcefield: {forcefield}, Water: {water_model}")
    print(f"{'='*60}")

    # Check if outputs already exist
    expected_outputs = [
        work_dir / "protein.gro",
        work_dir / "protein.top",
        work_dir / "posre.itp",
    ]
    if all(p.exists() for p in expected_outputs) and not force:
        print(f"✓ Outputs already exist. Skipping (use --force to re-run)")
        return work_dir

    # Step 1: Preprocess PDB for pH
    print(f"\nPreprocessing PDB for pH {ph}...")
    pdb_preprocessed = work_dir / "aSyn_preprocessed.pdb"
    preprocess_pdb_for_ph(pdb_src, pdb_preprocessed, ph_prot_cfg)
    print(f"  → {pdb_preprocessed}")

    # Step 2: Run pdb2gmx
    print(f"\nRunning pdb2gmx...")
    gro_out = work_dir / "protein.gro"
    top_out = work_dir / "protein.top"
    itp_out = work_dir / "posre.itp"

    run_pdb2gmx(
        pdb_preprocessed,
        gro_out,
        top_out,
        itp_out,
        gmx_path,
        forcefield,
        water_model,
    )

    print(f"  → {gro_out}")
    print(f"  → {top_out}")
    print(f"  → {itp_out}")

    # Step 3: Validation
    print(f"\nValidating outputs...")
    validate_protein_outputs(work_dir)

    print(f"\n✓ Protein topology at pH {ph} complete")
    return work_dir


def run_pdb2gmx(
    input_pdb: Path,
    gro_out: Path,
    top_out: Path,
    itp_out: Path,
    gmx_path: str,
    forcefield: str,
    water_model: str,
) -> None:
    """Run gmx pdb2gmx in non-interactive mode.

    Args:
        input_pdb: Input PDB file (pH-preprocessed)
        gro_out: Output .gro file
        top_out: Output .top file
        itp_out: Output position restraint .itp file
        gmx_path: Path to gmx executable
        forcefield: Force field name (e.g., 'amber19sb')
        water_model: Water model (e.g., 'opc')
    """
    cmd = [
        gmx_path,
        "pdb2gmx",
        "-f",
        str(input_pdb),
        "-o",
        str(gro_out),
        "-p",
        str(top_out),
        "-i",
        str(itp_out),
        "-ff",
        forcefield,
        "-water",
        water_model,
        "-ignh",
    ]

    work_dir = input_pdb.parent

    # Run with output capture for log file
    try:
        result = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        # Write logs
        log_file = work_dir / "pdb2gmx.log"
        with open(log_file, "w") as f:
            f.write(f"Command: {' '.join(cmd)}\n\n")
            f.write("=== STDOUT ===\n")
            f.write(e.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(e.stderr)

        raise RuntimeError(
            f"pdb2gmx failed with return code {e.returncode}\n"
            f"See {log_file} for details"
        )

    # Write log on success
    log_file = work_dir / "pdb2gmx.log"
    with open(log_file, "w") as f:
        f.write(f"Command: {' '.join(cmd)}\n\n")
        f.write("=== STDOUT ===\n")
        f.write(result.stdout)
        if result.stderr:
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)


def validate_protein_outputs(output_dir: Path) -> None:
    """Validate protein topology outputs.

    Args:
        output_dir: Directory containing protein.gro, protein.top, posre.itp
    """
    gro_file = output_dir / "protein.gro"
    top_file = output_dir / "protein.top"
    itp_file = output_dir / "posre.itp"

    # Check existence
    required_files = [gro_file, top_file, itp_file]
    missing = [f for f in required_files if not f.exists()]
    if missing:
        raise FileNotFoundError(f"Missing output files: {[str(f) for f in missing]}")

    # Validate .gro file
    with open(gro_file) as f:
        gro_lines = f.readlines()

    if len(gro_lines) < 3:
        raise ValueError(f".gro file too short: {len(gro_lines)} lines")

    # Count atoms (skip title, comment, atom count, and box vector lines)
    # Format: title, natoms, atoms..., box vector
    try:
        n_atoms = int(gro_lines[1].strip())
    except ValueError:
        raise ValueError(f"Cannot parse atom count from .gro file")

    if n_atoms < 900:
        raise ValueError(
            f".gro has only {n_atoms} atoms, expected >900 for aSyn + water + ions"
        )

    print(f"  ✓ .gro file: {n_atoms} atoms")

    # Validate .top file
    with open(top_file) as f:
        top_content = f.read()

    required_includes = ["amber19sb.ff/forcefield.itp"]
    for include in required_includes:
        if include not in top_content:
            raise ValueError(f".top missing include: {include}")

    water_model_line = None
    if "amber19sb.ff/opc.itp" in top_content:
        water_model_line = "opc"
    elif "amber19sb.ff/tip3p.itp" in top_content:
        water_model_line = "tip3p"
    elif "amber19sb.ff/tip4p.itp" in top_content:
        water_model_line = "tip4p"

    if water_model_line:
        print(f"  ✓ .top includes water model: {water_model_line}")
    else:
        raise ValueError(".top missing water model include")

    if "[ system ]" not in top_content:
        raise ValueError(".top missing [ system ] section")

    print(f"  ✓ .top file validated")

    # Validate posre.itp file
    with open(itp_file) as f:
        itp_content = f.read()

    if "[ position_restraints ]" not in itp_content:
        raise ValueError(".itp missing [ position_restraints ] section")

    print(f"  ✓ posre.itp file validated")

    # Check for pdb2gmx errors in log
    log_file = output_dir / "pdb2gmx.log"
    if log_file.exists():
        with open(log_file) as f:
            log_content = f.read()

        for line in log_content.split("\n"):
            if "CRITICAL" in line or "Fatal" in line:
                print(f"  ⚠ Warning: pdb2gmx warning detected: {line}")

    print(f"  ✓ All validations passed")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate GROMACS protein topologies via pdb2gmx"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/topology_config.yaml"),
        help="Path to topology config YAML",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if outputs already exist",
    )
    args = parser.parse_args()

    # Load config
    config = load_config(args.config)
    repo_root = Path(config["paths"]["repo_root"])
    os.chdir(repo_root)

    output_root = repo_root / config["paths"]["output_dir"]
    output_root.mkdir(parents=True, exist_ok=True)

    gmx_path = config["paths"]["gmx"]
    protein_cfg = config["protein"]
    conditions = config["conditions"]
    protonation_cfg = config["protonation"]

    # Convert relative path in config to absolute
    if not Path(protein_cfg["input_pdb"]).is_absolute():
        protein_cfg["input_pdb"] = repo_root / protein_cfg["input_pdb"]

    # Collect unique pH values
    unique_phs = set()
    for cond_cfg in conditions.values():
        unique_phs.add(cond_cfg["ph"])

    print(f"\nTopology Generation: Protein (pdb2gmx) Pipeline")
    print(f"Repository: {repo_root}")
    print(f"Output: {output_root}")
    print(f"Unique pH values: {sorted(unique_phs)}")
    print(f"GMX path: {gmx_path}")

    # Check GMX exists
    if not Path(gmx_path).exists():
        raise FileNotFoundError(f"GROMACS executable not found: {gmx_path}")

    # Process each unique pH
    for ph in sorted(unique_phs):
        ph_str = str(ph)
        if ph_str not in protonation_cfg:
            raise ValueError(
                f"pH {ph} in conditions but not in protonation config"
            )

        ph_prot_cfg = protonation_cfg[ph_str]

        process_ph(
            ph,
            protein_cfg,
            ph_prot_cfg,
            output_root,
            gmx_path,
            force=args.force,
        )

    print(f"\n{'='*60}")
    print(f"✓ All protein topologies generated successfully")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
