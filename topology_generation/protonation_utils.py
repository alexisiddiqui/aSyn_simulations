"""pdb2pqr-based protonation utilities for the topology generation pipeline.

Replaces the naive string-based `preprocess_pdb_for_ph` with a proper
PropKa-driven protonation step that correctly handles all titratable
residues (HIS, ASP, GLU, LYS, CYS) at arbitrary pH values.

Workflow:
    input.pdb
        -> run_pdb2pqr()       (PropKa + AMBER protonation; outputs clean PDB)
        -> apply_manual_overrides()  (force per-residue overrides from config)
        -> pdb2gmx -ignh       (strips H, re-adds from FF, generates topology)
"""

import subprocess
import sys
from pathlib import Path
from typing import Optional


def run_pdb2pqr(
    input_pdb: Path,
    output_pdb: Path,
    ph: float,
    work_dir: Path,
) -> Path:
    """Run pdb2pqr with PropKa to assign protonation states at a given pH.

    Uses the AMBER force field naming convention so residue names (HIE, HID,
    HIP, GLH, ASH, etc.) are directly compatible with ``gmx pdb2gmx -ff amber*``.

    A PQR file is also written alongside the PDB as
    ``<output_pdb.stem>.pqr`` — it is kept for reference but not used
    downstream (``pdb2gmx`` uses the PDB).

    Args:
        input_pdb: Raw input PDB (no explicit hydrogens required; pdb2pqr adds
                   them, and pdb2gmx will strip + re-add from the FF).
        output_pdb: Destination for the protonated, AMBER-named PDB file.
        ph: Target pH for PropKa titration (e.g. 7.4, 7.2, 4.9).
        work_dir: Working directory; log is written to ``pdb2pqr.log`` here.

    Returns:
        Path to the protonated PDB file (same as ``output_pdb``).

    Raises:
        RuntimeError: If pdb2pqr exits with a non-zero return code.
        FileNotFoundError: If the output PDB was not created despite a zero
                           return code (should not happen in practice).
    """
    input_pdb = Path(input_pdb)
    output_pdb = Path(output_pdb)
    work_dir = Path(work_dir)

    output_pqr = output_pdb.with_suffix(".pqr")

    # Resolve the pdb2pqr30 executable relative to the current Python interpreter
    # so that the script works whether or not the venv is activated in the shell.
    _bin_dir = Path(sys.executable).parent
    _pdb2pqr_bin = _bin_dir / "pdb2pqr30"
    pdb2pqr_exe = str(_pdb2pqr_bin) if _pdb2pqr_bin.exists() else "pdb2pqr30"

    cmd = [
        pdb2pqr_exe,
        "--ff", "AMBER",
        "--with-ph", str(ph),
        "--titration-state-method", "propka",
        "--pdb-output", str(output_pdb),
        "--nodebump",          # don't move heavy atoms to fix steric clashes
        "--noopt",             # don't optimise H-bond network (pdb2gmx does this)
        str(input_pdb),
        str(output_pqr),
    ]

    log_file = work_dir / "pdb2pqr.log"

    try:
        result = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        with open(log_file, "w") as f:
            f.write(f"Command: {' '.join(cmd)}\n\n")
            f.write("=== STDOUT ===\n")
            f.write(e.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(e.stderr)
        raise RuntimeError(
            f"pdb2pqr failed (return code {e.returncode}) at pH {ph}.\n"
            f"See {log_file} for details."
        )

    # Write log on success
    with open(log_file, "w") as f:
        f.write(f"Command: {' '.join(cmd)}\n\n")
        f.write("=== STDOUT ===\n")
        f.write(result.stdout)
        if result.stderr:
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)

    if not output_pdb.exists():
        raise FileNotFoundError(
            f"pdb2pqr reported success but output PDB not found: {output_pdb}\n"
            f"See {log_file} for details."
        )

    return output_pdb


def apply_manual_overrides(pdb_path: Path, overrides: Optional[dict]) -> None:
    """Apply per-residue protonation overrides to a pdb2pqr-generated PDB.

    Called after ``run_pdb2pqr`` to enforce any entries in the ``overrides``
    dict from ``topology_config.yaml`` (e.g. ``"A:50": "HIP"``).  Edits the
    file in-place.

    Args:
        pdb_path: Path to the protonated PDB (output of ``run_pdb2pqr``).
        overrides: Dict mapping ``"<chain>:<resnum>"`` → new residue name.
                   Pass ``None`` or ``{}`` to skip (no-op).
    """
    if not overrides:
        return

    pdb_path = Path(pdb_path)

    with open(pdb_path) as f:
        lines = f.readlines()

    output_lines = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            # PDB fixed-width columns (1-indexed as per spec):
            #   18-20  residue name   (0-indexed 17:20)
            #   22     chain ID       (0-indexed 21)
            #   23-26  residue seq    (0-indexed 22:26)
            chain = line[21].strip()
            resnum = line[22:26].strip()
            key = f"{chain}:{resnum}"

            if key in overrides:
                new_resname = overrides[key]
                # Right-justify in a 3-char field (GROMACS convention)
                line = line[:17] + new_resname.rjust(3) + line[20:]

        output_lines.append(line)

    with open(pdb_path, "w") as f:
        f.writelines(output_lines)


def summarise_protonation(pdb_path: Path) -> dict[str, list[str]]:
    """Return a dict of titratable residue names found in the PDB.

    Useful for logging and validation.  Returns a mapping of residue name
    to a list of ``"<chain>:<resnum>"`` occurrences.

    Args:
        pdb_path: Protonated PDB file to inspect.
    """
    from collections import defaultdict

    pdb_path = Path(pdb_path)
    titratable = {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX", "CYM"}
    found: dict[str, list[str]] = defaultdict(list)

    seen: set[str] = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = line[22:26].strip()
            key = f"{chain}:{resnum}"
            if resname in titratable and key not in seen:
                found[resname].append(key)
                seen.add(key)

    return dict(found)
