"""Shared utilities for GROMACS system building pipeline."""

import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import yaml


# ---------------------------------------------------------------------------
# Count calculations
# ---------------------------------------------------------------------------

def compute_box_volume(size_nm: float, shape: str) -> float:
    """Compute box volume in nm³ from box vector length and shape.

    Args:
        size_nm: Box vector length L in nm
        shape: 'cubic' or 'dodecahedron'

    Returns:
        Volume in nm³
    """
    if shape == "cubic":
        return size_nm ** 3
    elif shape == "dodecahedron":
        return size_nm ** 3 / math.sqrt(2)
    else:
        raise ValueError(f"Unknown box shape: {shape!r}. Use 'cubic' or 'dodecahedron'")


def compute_molecule_count(conc_mm: float, volume_nm3: float, min_count: int = 0) -> int:
    """Compute number of molecules for a given concentration and volume.

    Formula: n = C_mM × V_nm³ × 6.022e-4
    Derivation: n = (C_mM×1e-3 mol/L) × (V_nm³×1e-24 L) × 6.022e23 /mol

    Args:
        conc_mm: Concentration in mM
        volume_nm3: Box volume in nm³
        min_count: Minimum number to return (0 = skip if rounds to 0)

    Returns:
        Number of molecules (>= min_count)
    """
    n = round(conc_mm * volume_nm3 * 6.022e-4)
    return max(n, min_count)


def compute_buffer_counts(
    buffer_species: dict,
    concentration_mm: float,
    volume_nm3: float,
) -> dict:
    """Compute molecule counts for buffer species from stoichiometric ratios.

    Args:
        buffer_species: Dict of {resname: ratio} e.g. {'TR1': 1.0, 'TR0': 0.22}
        concentration_mm: Total buffer concentration in mM
        volume_nm3: Box volume in nm³

    Returns:
        Dict of {resname: count}, preserving config order
    """
    if not buffer_species:
        return {}

    n_total = compute_molecule_count(concentration_mm, volume_nm3)
    species = list(buffer_species.items())

    if len(species) == 1:
        return {species[0][0]: n_total}

    ratio_sum = sum(r for _, r in species)
    counts = {}
    remaining = n_total

    for resname, ratio in species[:-1]:
        n = round(n_total * ratio / ratio_sum)
        counts[resname] = n
        remaining -= n

    # Last species gets the remainder to ensure counts sum exactly to n_total
    counts[species[-1][0]] = max(remaining, 0)
    return counts


def compute_ion_counts(
    ions_cfg: dict,
    volume_nm3: float,
) -> list[tuple[str, int, int]]:
    """Compute ion molecule counts, sorted largest-first.

    Args:
        ions_cfg: Dict from system_config conditions[x].ions
                  e.g. {'NA': {conc_mm: 143, charge: 1}, ...}
        volume_nm3: Box volume in nm³

    Returns:
        List of (ion_name, count, charge) sorted by count descending, skipping zeros
    """
    results = []
    for ion_name, cfg in ions_cfg.items():
        min_count = cfg.get("min_count", 0)
        n = compute_molecule_count(cfg["conc_mm"], volume_nm3, min_count=min_count)
        if n > 0:
            results.append((ion_name, n, cfg["charge"]))

    # Sort largest count first for genion efficiency
    results.sort(key=lambda x: x[1], reverse=True)
    return results


# ---------------------------------------------------------------------------
# Topology assembly
# ---------------------------------------------------------------------------

def _extract_gaff2_atomtypes(top_paths: list) -> str:
    """Extract and deduplicate [ atomtypes ] lines from acpype .top files.

    Args:
        top_paths: List of Path objects pointing to {resname}_GMX.top files

    Returns:
        String with deduplicated atomtype lines (without the [ atomtypes ] header)
    """
    seen = {}  # name -> full line (first occurrence wins)

    for top_path in top_paths:
        with open(top_path) as f:
            content = f.read()

        start = content.find("[ atomtypes ]")
        if start == -1:
            continue

        in_section = False
        for line in content[start:].splitlines():
            if line.strip() == "[ atomtypes ]":
                in_section = True
                continue
            if in_section:
                stripped = line.strip()
                if stripped.startswith("["):
                    break
                if not stripped or stripped.startswith(";"):
                    continue
                name = stripped.split()[0]
                if name not in seen:
                    seen[name] = line

    return "\n".join(seen.values())


def _clean_ligand_itp(itp_path: Path, work_dir: Path) -> Path:
    """Copy ligand ITP to work_dir, stripping [ system ] and [ molecules ] sections.

    acpype includes these sections in .itp files; they must be removed before
    including the ITP in a larger topology.

    Args:
        itp_path: Source {resname}_GMX.itp
        work_dir: Destination directory

    Returns:
        Path to cleaned copy in work_dir
    """
    with open(itp_path) as f:
        content = f.read()

    # Truncate at [ system ] (preceded by newline)
    for marker in ["\n[ system ]", "\n [ system ]"]:
        idx = content.find(marker)
        if idx != -1:
            content = content[:idx] + "\n"
            break

    clean_path = work_dir / itp_path.name
    with open(clean_path, "w") as f:
        f.write(content)

    return clean_path


def assemble_topology(
    protein_top: Path,
    ligand_info: list,
    buffer_counts: dict,
    work_dir: Path,
) -> Path:
    """Assemble system.top from protein.top + ligand ITPs.

    Modifications made to protein.top copy:
    1. Insert combined GAFF2 [ atomtypes ] ITP after forcefield include
    2. Insert cleaned ligand molecule ITPs before [ system ]
    3. Append buffer molecule counts to [ molecules ]

    Args:
        protein_top: Path to protein.top (from pdb2gmx)
        ligand_info: List of (resname, itp_path, top_path) for each ligand in condition
        buffer_counts: Dict of {resname: count} to add to [ molecules ]
        work_dir: Output directory for system.top and helper files

    Returns:
        Path to assembled system.top
    """
    system_top = work_dir / "system.top"

    # Build combined GAFF2 atomtypes ITP
    atomtypes_itp = None
    if ligand_info:
        top_paths = [top for _, _, top in ligand_info]
        atomtypes_lines = _extract_gaff2_atomtypes(top_paths)
        if atomtypes_lines:
            atomtypes_itp = work_dir / "gaff2_atomtypes.itp"
            with open(atomtypes_itp, "w") as f:
                f.write("; Combined GAFF2 atomtypes — generated by build_system.py\n")
                f.write("; nbfunc=1 comb-rule=2 (AMBER/GAFF2 convention)\n\n")
                f.write("[ atomtypes ]\n")
                f.write(";name   bond_type     mass     charge   ptype   sigma         epsilon\n")
                f.write(atomtypes_lines + "\n")

    # Create cleaned ligand ITP copies (strip [ system ] / [ molecules ])
    clean_itps = []  # [(resname, Path)]
    for resname, itp_path, _ in ligand_info:
        clean_path = _clean_ligand_itp(itp_path, work_dir)
        clean_itps.append((resname, clean_path))

    # Read protein.top
    with open(protein_top) as f:
        content = f.read()

    # 1. Insert atomtypes include after forcefield include
    ff_include = '#include "amber19sb.ff/forcefield.itp"'
    if atomtypes_itp and ff_include in content:
        insert = f'\n#include "{atomtypes_itp.absolute()}"\n; end GAFF2 atomtypes'
        content = content.replace(ff_include, ff_include + insert, 1)

    # 2. Insert ligand molecule ITPs before [ system ]
    if clean_itps:
        ligand_block = "\n; Ligand molecule topologies\n"
        for resname, clean_path in clean_itps:
            ligand_block += f'#include "{clean_path.absolute()}"\n'
        ligand_block += "\n"
        content = content.replace("[ system ]", ligand_block + "[ system ]", 1)

    # 3. Append buffer molecule counts after Protein_chain_A in [ molecules ]
    non_zero = {r: n for r, n in buffer_counts.items() if n > 0}
    if non_zero:
        buffer_lines = "\n".join(
            f"{resname:<20} {count}" for resname, count in non_zero.items()
        )
        # Only replace within the [ molecules ] section (not [ moleculetype ])
        mol_section_idx = content.rfind("[ molecules ]")
        if mol_section_idx == -1:
            raise ValueError(f"No [ molecules ] section found in {protein_top}")

        mol_section = content[mol_section_idx:]
        pattern = re.compile(r"(Protein_chain_A\s+\d+)", re.MULTILINE)
        if pattern.search(mol_section):
            mol_section = pattern.sub(r"\1\n" + buffer_lines, mol_section, count=1)
            content = content[:mol_section_idx] + mol_section
        else:
            raise ValueError(
                f"Could not find 'Protein_chain_A' in [ molecules ] section of {protein_top}"
            )

    with open(system_top, "w") as f:
        f.write(content)

    return system_top


# ---------------------------------------------------------------------------
# GROMACS command wrappers
# ---------------------------------------------------------------------------

def run_gmx(
    gmx_path: str,
    args: list,
    cwd: Path,
    label: str,
    stdin: Optional[str] = None,
) -> None:
    """Run a GROMACS command, writing stdout/stderr to {label}.log.

    Args:
        gmx_path: Path to gmx executable
        args: Arguments after 'gmx' (e.g. ['editconf', '-f', 'in.gro', ...])
        cwd: Working directory
        label: Log file prefix ({label}.log written to cwd)
        stdin: Optional string piped to stdin (e.g. 'SOL\\n' for genion)

    Raises:
        RuntimeError: If command exits non-zero
    """
    cmd = [gmx_path] + args
    log_file = Path(cwd) / f"{label}.log"

    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            check=True,
            input=stdin,
        )
    except subprocess.CalledProcessError as e:
        with open(log_file, "w") as f:
            f.write(f"Command: {' '.join(cmd)}\n\n")
            f.write("=== STDOUT ===\n")
            f.write(e.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(e.stderr)
        raise RuntimeError(
            f"GROMACS command failed: {' '.join(cmd)}\n"
            f"Return code: {e.returncode}\n"
            f"See {log_file} for details"
        )

    with open(log_file, "w") as f:
        f.write(f"Command: {' '.join(cmd)}\n\n")
        f.write("=== STDOUT ===\n")
        f.write(result.stdout)
        if result.stderr:
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)


def run_editconf(
    gmx_path: str,
    input_gro: Path,
    output_gro: Path,
    box_size: float,
    box_shape: str,
    work_dir: Path,
) -> None:
    """Center protein and define simulation box.

    Args:
        gmx_path: Path to gmx executable
        input_gro: Input structure (protein.gro)
        output_gro: Output structure with box defined
        box_size: Box vector length in nm
        box_shape: 'cubic' or 'dodecahedron'
        work_dir: Working directory for log
    """
    run_gmx(
        gmx_path,
        [
            "editconf",
            "-f", str(input_gro),
            "-o", str(output_gro),
            "-c",
            "-box", str(box_size), str(box_size), str(box_size),
            "-bt", box_shape,
            "-princ"
        ],
        work_dir,
        "editconf",
    )


def run_insert_molecules(
    gmx_path: str,
    input_gro: Path,
    molecule_gro: Path,
    n_mols: int,
    output_gro: Path,
    work_dir: Path,
    label: str,
) -> int:
    """Insert molecules into a box.

    Args:
        gmx_path: Path to gmx executable
        input_gro: Current system structure
        molecule_gro: Single-molecule structure to insert
        n_mols: Number of molecules to insert
        output_gro: Output structure
        work_dir: Working directory for log
        label: Log file prefix

    Returns:
        Number of molecules actually inserted (parsed from log)
    """
    run_gmx(
        gmx_path,
        [
            "insert-molecules",
            "-f", str(input_gro),
            "-ci", str(molecule_gro),
            "-nmol", str(n_mols),
            "-o", str(output_gro),
        ],
        work_dir,
        label,
    )

    # Parse inserted count from log
    log_file = work_dir / f"{label}.log"
    inserted = n_mols  # default assumption
    if log_file.exists():
        content = log_file.read_text()
        for line in content.splitlines():
            if "Added" in line and "molecules" in line:
                try:
                    inserted = int(line.split()[1])
                except (IndexError, ValueError):
                    pass

    if inserted != n_mols:
        print(
            f"  WARNING: Requested {n_mols} molecules but inserted {inserted}. "
            f"Box may be too crowded."
        )

    return inserted


def run_solvate(
    gmx_path: str,
    input_gro: Path,
    top_path: Path,
    output_gro: Path,
    work_dir: Path,
    water_model: str,
) -> None:
    """Solvate the system with a model-aware solvent box.

    Uses spc216.gro for 3-site models (opc3, tip3p, etc.)
    and tip4p.gro for 4-site models (opc, tip4p).

    Args:
        gmx_path: Path to gmx executable
        input_gro: Current system structure (protein + ligands)
        top_path: Topology file (updated in-place with SOL count)
        output_gro: Solvated output structure
        work_dir: Working directory for log
        water_model: Name of water model (e.g., 'opc3', 'tip3p', 'opc')
    """
    # Select solvent box based on sites
    three_site_models = ["opc3", "tip3p", "spce", "spc", "tip3p-fb"]
    four_site_models = ["opc", "tip4p", "tip4p-fb", "tip4p-ew"]

    if water_model.lower() in three_site_models:
        solvent_box = "spc216.gro"
    elif water_model.lower() in four_site_models:
        solvent_box = "tip4p.gro"
    else:
        # Fallback to spc216 for unknown models, or we could raise error
        print(f"  WARNING: Unknown water model {water_model!r}, defaulting to spc216.gro")
        solvent_box = "spc216.gro"

    run_gmx(
        gmx_path,
        [
            "solvate",
            "-cp", str(input_gro),
            "-cs", solvent_box,
            "-o", str(output_gro),
            "-p", str(top_path),
        ],
        work_dir,
        "solvate",
    )


def run_grompp_for_genion(
    gmx_path: str,
    ions_mdp: Path,
    input_gro: Path,
    top_path: Path,
    output_tpr: Path,
    work_dir: Path,
    label: str,
) -> None:
    """Run grompp to generate .tpr for genion.

    Args:
        gmx_path: Path to gmx executable
        ions_mdp: Path to ions.mdp
        input_gro: Current system structure
        top_path: Topology file
        output_tpr: Output .tpr file
        work_dir: Working directory for log
        label: Log file prefix
    """
    run_gmx(
        gmx_path,
        [
            "grompp",
            "-f", str(ions_mdp),
            "-c", str(input_gro),
            "-p", str(top_path),
            "-o", str(output_tpr),
            "-maxwarn", "2",
        ],
        work_dir,
        label,
    )


def run_genion(
    gmx_path: str,
    tpr_path: Path,
    top_path: Path,
    output_gro: Path,
    work_dir: Path,
    label: str,
    n_pos: int = 0,
    pname: str = "NA",
    pcharge: int = 1,
    n_neg: int = 0,
    nname: str = "CL",
    ncharge: int = -1,
    neutral: bool = False,
) -> None:
    """Run genion to add monoatomic ions.

    Args:
        gmx_path: Path to gmx executable
        tpr_path: Input .tpr from grompp
        top_path: Topology file (updated in-place)
        output_gro: Output structure with ions
        work_dir: Working directory for log
        label: Log file prefix
        n_pos: Number of positive ions to add
        pname: Positive ion molecule name (e.g. 'NA', 'K', 'CA')
        pcharge: Charge of positive ion
        n_neg: Number of negative ions to add
        nname: Negative ion molecule name (e.g. 'CL')
        ncharge: Charge of negative ion (negative int)
        neutral: If True, add extra ions to neutralize system charge
    """
    args = [
        "genion",
        "-s", str(tpr_path),
        "-o", str(output_gro),
        "-p", str(top_path),
        "-pname", pname,
        "-pq", str(pcharge),
        "-nname", nname,
        "-nq", str(ncharge),
    ]

    if n_pos > 0:
        args += ["-np", str(n_pos)]
    if n_neg > 0:
        args += ["-nn", str(n_neg)]
    if neutral:
        args.append("-neutral")

    # genion prompts for solvent group interactively; pipe "SOL" to stdin
    run_gmx(gmx_path, args, work_dir, label, stdin="SOL\n")
