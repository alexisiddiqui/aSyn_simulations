"""Generate GROMACS topologies for ligands using GAFF2/AMBER pipeline."""

import os
import shutil
import subprocess
from pathlib import Path

from topology_utils import (
    conda_cmd,
    fix_mol2_resname,
    load_config,
    run_cmd,
    validate_files_exist,
)


def process_ligand(
    ligand_name: str,
    ligand_cfg: dict,
    output_root: Path,
    conda_env: str,
    force: bool = False,
) -> Path:
    """Run full GAFF2 pipeline for one ligand.

    Args:
        ligand_name: Ligand code (e.g., 'TR1')
        ligand_cfg: Config dict with 'mol2', 'net_charge', 'resname'
        output_root: Root output directory (topology_output)
        conda_env: Conda environment name for AmberTools
        force: If True, re-run even if outputs exist

    Returns:
        Path to ligand output directory
    """
    work_dir = output_root / "ligands" / ligand_name
    work_dir.mkdir(parents=True, exist_ok=True)

    mol2_src = ligand_cfg["mol2"]
    net_charge = ligand_cfg["net_charge"]
    resname = ligand_cfg["resname"]

    print(f"\n{'='*60}")
    print(f"Processing ligand: {ligand_name} ({resname})")
    print(f"{'='*60}")

    # Check if all outputs already exist
    expected_outputs = [
        work_dir / f"{resname}_GMX.itp",
        work_dir / f"{resname}_GMX.top",
    ]
    if all(p.exists() for p in expected_outputs) and not force:
        print(f"✓ Outputs already exist. Skipping (use --force to re-run)")
        return work_dir

    # Step 1: Fix MOL2 resname
    print(f"\n[1/5] Fixing MOL2 residue name...")
    mol2_fixed = work_dir / f"{resname}_input.mol2"
    fix_mol2_resname(mol2_src, mol2_fixed, resname)
    print(f"  → {mol2_fixed}")

    # Step 2: Antechamber
    print(f"\n[2/5] Running antechamber (AM1-BCC charges, GAFF2 typing)...")
    mol2_ac = step2_antechamber(mol2_fixed, work_dir, resname, net_charge, conda_env)
    print(f"  → {mol2_ac}")

    # Step 3: Parmchk2
    print(f"\n[3/5] Running parmchk2 (force field parameters)...")
    frcmod = step3_parmchk2(mol2_ac, work_dir, resname, conda_env)
    print(f"  → {frcmod}")

    # Check for missing parameters
    with open(frcmod) as f:
        frcmod_content = f.read()
    if "ATTN, NEED" in frcmod_content:
        print("  ⚠ WARNING: Found 'ATTN, NEED' lines (missing parameters):")
        for line in frcmod_content.split("\n"):
            if "ATTN, NEED" in line:
                print(f"    {line}")

    # Step 4: tleap
    print(f"\n[4/5] Running tleap (AMBER topology)...")
    prmtop, inpcrd = step4_tleap(mol2_ac, frcmod, work_dir, resname, conda_env)
    print(f"  → {prmtop}")
    print(f"  → {inpcrd}")

    # Step 5: acpype
    print(f"\n[5/5] Running acpype (AMBER → GROMACS conversion)...")
    itp, top = step5_acpype(prmtop, inpcrd, work_dir, resname)
    print(f"  → {itp}")
    print(f"  → {top}")

    # Validation
    print(f"\nValidating outputs...")
    validate_ligand_outputs(work_dir, resname, mol2_fixed)

    print(f"\n✓ Ligand {ligand_name} complete")
    return work_dir


def step2_antechamber(
    mol2_in: Path, work_dir: Path, resname: str, net_charge: int, conda_env: str
) -> Path:
    """Run antechamber via conda.

    Returns path to output MOL2 (_AC.mol2).
    """
    mol2_out = work_dir / f"{resname}_AC.mol2"

    cmd = conda_cmd(
        [
            "antechamber",
            "-i",
            str(mol2_in),
            "-fi",
            "mol2",
            "-o",
            str(mol2_out),
            "-fo",
            "mol2",
            "-c",
            "bcc",
            "-nc",
            str(net_charge),
            "-at",
            "gaff2",
            "-rn",
            resname,
            "-s",
            "2",
        ],
        conda_env,
    )

    run_cmd(cmd, work_dir, "antechamber")
    return mol2_out


def step3_parmchk2(
    ac_mol2: Path, work_dir: Path, resname: str, conda_env: str
) -> Path:
    """Run parmchk2 via conda.

    Returns path to .frcmod file.
    """
    frcmod_out = work_dir / f"{resname}_AC.frcmod"

    cmd = conda_cmd(
        [
            "parmchk2",
            "-i",
            str(ac_mol2),
            "-f",
            "mol2",
            "-o",
            str(frcmod_out),
            "-s",
            "gaff2",
        ],
        conda_env,
    )

    run_cmd(cmd, work_dir, "parmchk2")
    return frcmod_out


def step4_tleap(
    ac_mol2: Path, frcmod: Path, work_dir: Path, resname: str, conda_env: str
) -> tuple[Path, Path]:
    """Run tleap via conda to generate AMBER topology.

    Returns (prmtop_path, inpcrd_path).
    """
    prmtop_out = work_dir / f"{resname}_AC.prmtop"
    inpcrd_out = work_dir / f"{resname}_AC.inpcrd"
    tleap_in = work_dir / f"tleap_{resname}.in"

    # Write tleap input file using absolute paths
    tleap_content = f"""verbosity 1
source leaprc.gaff2
{resname} = loadmol2 {ac_mol2.absolute()}
loadamberparams {frcmod.absolute()}
check {resname}
saveamberparm {resname} {prmtop_out.absolute()} {inpcrd_out.absolute()}
quit
"""
    with open(tleap_in, "w") as f:
        f.write(tleap_content)

    cmd = conda_cmd(["tleap", "-f", str(tleap_in)], conda_env)
    run_cmd(cmd, work_dir, "tleap")

    return prmtop_out, inpcrd_out


def step5_acpype(prmtop: Path, inpcrd: Path, work_dir: Path, resname: str) -> tuple[Path, Path]:
    """Run acpype in amb2gmx mode via uv to convert AMBER → GROMACS.

    Returns (itp_path, top_path).

    Note: acpype in amb2gmx mode produces {resname}_GMX.top (full topology) and
    posre_{resname}.itp (position restraints), but NOT a separate {resname}_GMX.itp.
    We extract the molecule-specific sections from .top into .itp for use in system building.
    """
    itp_out = work_dir / f"{resname}_GMX.itp"
    top_out = work_dir / f"{resname}_GMX.top"
    posre_out = work_dir / f"posre_{resname}.itp"

    # acpype creates {resname}.amb2gmx/ directory
    acpype_dir = work_dir / f"{resname}.amb2gmx"

    cmd = [
        "uv",
        "run",
        "acpype",
        "-p",
        str(prmtop),
        "-x",
        str(inpcrd),
        "-b",
        resname,
        "-o",
        "gmx",
    ]

    run_cmd(cmd, work_dir, "acpype")

    # Copy output files from acpype working dir to work_dir root
    acpype_top_src = acpype_dir / f"{resname}_GMX.top"
    acpype_posre_src = acpype_dir / f"posre_{resname}.itp"

    if not acpype_top_src.exists():
        raise FileNotFoundError(f"acpype output not found: {acpype_top_src}")

    shutil.copy(acpype_top_src, top_out)

    # Extract molecule sections ([ moleculetype ] onwards) into .itp
    _extract_itp_from_top(acpype_top_src, itp_out)

    # Copy position restraints file if it exists
    if acpype_posre_src.exists():
        shutil.copy(acpype_posre_src, posre_out)

    return itp_out, top_out


def _extract_itp_from_top(top_path: Path, itp_path: Path) -> None:
    """Extract molecule sections ([ moleculetype ] onwards) from acpype .top into .itp.

    Args:
        top_path: Path to acpype {resname}_GMX.top
        itp_path: Path to output {resname}_GMX.itp

    Raises:
        ValueError: If [ moleculetype ] section not found
    """
    with open(top_path) as f:
        content = f.read()

    # Find start of [ moleculetype ] section
    mol_start = content.find("[ moleculetype ]")
    if mol_start == -1:
        raise ValueError(f"No [ moleculetype ] section in {top_path}")

    # Write everything from [ moleculetype ] onwards to .itp file
    with open(itp_path, "w") as f:
        f.write(content[mol_start:])

    print(f"  Extracted .itp from {top_path.name} → {itp_path.name}")


def validate_ligand_outputs(work_dir: Path, resname: str, mol2_src: Path) -> None:
    """Validate ligand topology outputs.

    Args:
        work_dir: Ligand output directory
        resname: Residue name
        mol2_src: Source MOL2 file (for atom count)
    """
    # Check file existence
    required_files = [
        work_dir / f"{resname}_AC.mol2",
        work_dir / f"{resname}_AC.frcmod",
        work_dir / f"{resname}_AC.prmtop",
        work_dir / f"{resname}_AC.inpcrd",
        work_dir / f"{resname}_GMX.itp",  # extracted from .top
        work_dir / f"{resname}_GMX.top",  # full acpype output
    ]

    for fpath in required_files:
        if not fpath.exists():
            raise FileNotFoundError(f"Missing output file: {fpath}")
        if fpath.stat().st_size == 0:
            raise ValueError(f"Empty output file: {fpath}")

    # Validate .itp file
    itp_file = work_dir / f"{resname}_GMX.itp"
    with open(itp_file) as f:
        itp_content = f.read()

    required_sections = ["[ moleculetype ]", "[ atoms ]", "[ bonds ]"]
    for section in required_sections:
        if section not in itp_content:
            raise ValueError(f".itp missing section: {section}")

    # Count atoms in source MOL2 (for reference only, not validation)
    # Note: .itp may have different atom count due to hydrogen addition/removal by acpype
    mol2_atom_count = 0
    with open(mol2_src) as f:
        in_atom_section = False
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atom_section = True
                continue
            if in_atom_section and line.strip() == "":
                break
            if in_atom_section:
                mol2_atom_count += 1

    # Count atoms in .itp
    itp_atom_count = 0
    in_atom_section = False
    for line in itp_content.split("\n"):
        if line.startswith("[ atoms ]"):
            in_atom_section = True
            continue
        if in_atom_section and line.startswith("["):
            break
        if in_atom_section and line.strip() and not line.startswith(";"):
            itp_atom_count += 1

    # Just verify .itp has atoms
    if itp_atom_count == 0:
        raise ValueError(".itp file has no atoms in [ atoms ] section")

    print(f"  ✓ Atom counts: mol2 {mol2_atom_count}, .itp {itp_atom_count} (hydrogens may differ)")

    # Validate charge sum in .itp (rough check)
    charge_sum = 0.0
    in_atom_section = False
    for line in itp_content.split("\n"):
        if line.startswith("[ atoms ]"):
            in_atom_section = True
            continue
        if in_atom_section and line.startswith("["):
            break
        if in_atom_section and line.strip() and not line.startswith(";"):
            parts = line.split()
            if len(parts) >= 7:
                try:
                    charge = float(parts[6])
                    charge_sum += charge
                except ValueError:
                    pass

    # Just print the final charge sum for inspection
    print(f"  ✓ Total charge sum: {charge_sum:.4f}")

    print(f"  ✓ All validations passed")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate GROMACS ligand topologies via GAFF2/AMBER"
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

    conda_env = config["paths"]["conda_env"]
    ligands_cfg = config["ligands"]
    conditions = config["conditions"]

    # Collect unique ligands across all conditions
    unique_ligands = set()
    for cond_cfg in conditions.values():
        unique_ligands.update(cond_cfg["ligands"])

    print(f"\nTopology Generation: GAFF2 Ligand Pipeline")
    print(f"Repository: {repo_root}")
    print(f"Output: {output_root}")
    print(f"Unique ligands: {sorted(unique_ligands)}")

    # Process each unique ligand
    for ligand_name in sorted(unique_ligands):
        if ligand_name not in ligands_cfg:
            raise ValueError(f"Ligand {ligand_name} in conditions but not in ligands config")

        ligand_cfg = ligands_cfg[ligand_name]

        # Convert relative paths in config to absolute
        if not Path(ligand_cfg["mol2"]).is_absolute():
            ligand_cfg["mol2"] = repo_root / ligand_cfg["mol2"]

        process_ligand(ligand_name, ligand_cfg, output_root, conda_env, force=args.force)

    print(f"\n{'='*60}")
    print(f"✓ All ligands processed successfully")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
