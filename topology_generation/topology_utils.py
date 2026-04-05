"""Shared utilities for topology generation pipeline."""

import subprocess
import sys
from pathlib import Path
from typing import Optional

import yaml


def load_config(config_path: Path) -> dict:
    """Load and validate YAML config.

    Args:
        config_path: Path to topology_config.yaml

    Returns:
        Parsed config dict

    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config is malformed
        KeyError: If required keys are missing
    """
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Validate required top-level keys
    required_keys = {"paths", "protein", "conditions", "ligands", "protonation"}
    missing = required_keys - set(config.keys())
    if missing:
        raise KeyError(f"Missing required config keys: {missing}")

    # Validate paths
    required_paths = {"repo_root", "gmx", "output_dir", "conda_env"}
    missing = required_paths - set(config["paths"].keys())
    if missing:
        raise KeyError(f"Missing required path keys: {missing}")

    return config


def run_cmd(cmd: list[str], cwd: Path, label: str, env: Optional[dict] = None) -> None:
    """Run a shell command and log output.

    Args:
        cmd: Command as list (e.g., ['echo', 'hello'])
        cwd: Working directory for command
        label: Label for log file (output goes to {label}.log in cwd)
        env: Optional environment dict (merged with os.environ)

    Raises:
        RuntimeError: If command fails (non-zero return code)
    """
    cwd = Path(cwd)
    log_file = cwd / f"{label}.log"

    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            check=True,
            env=env,
        )
    except subprocess.CalledProcessError as e:
        # Write logs even on failure for debugging
        with open(log_file, "w") as f:
            f.write(f"Command: {' '.join(cmd)}\n\n")
            f.write("=== STDOUT ===\n")
            f.write(e.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(e.stderr)

        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"Return code: {e.returncode}\n"
            f"See {log_file} for details"
        )

    # Write log on success
    with open(log_file, "w") as f:
        f.write(f"Command: {' '.join(cmd)}\n\n")
        f.write("=== STDOUT ===\n")
        f.write(result.stdout)
        if result.stderr:
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)


def fix_mol2_resname(src: Path, dst: Path, new_resname: str) -> None:
    """Rewrite MOL2 residue name in @<TRIPOS>MOLECULE and @<TRIPOS>ATOM sections.

    Args:
        src: Input MOL2 file
        dst: Output MOL2 file
        new_resname: New 3-letter residue name (left-justified, padded to 3 chars)
    """
    src = Path(src)
    dst = Path(dst)

    with open(src) as f:
        lines = f.readlines()

    new_resname_padded = new_resname.ljust(3)
    in_atom_section = False
    output_lines = []

    for i, line in enumerate(lines):
        # Detect @<TRIPOS>MOLECULE section (usually line 0-1)
        if line.startswith("@<TRIPOS>MOLECULE"):
            in_atom_section = False
            output_lines.append(line)
            # Next line is the molecule name; rewrite it with the new resname
            if i + 1 < len(lines):
                # Skip one more line to get to the name, then rewrite it
                i_next = i + 1
                # The molecule name line is typically just resname + some suffix
                # Simplest: just replace with resname
                continue
        elif i > 0 and lines[i-1].startswith("@<TRIPOS>MOLECULE"):
            # This is the molecule name line; rewrite with new resname
            output_lines.append(f"{new_resname_padded}\n")
        elif line.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            output_lines.append(line)
        elif in_atom_section and line.strip() == "":
            in_atom_section = False
            output_lines.append(line)
        elif in_atom_section and line.strip():
            # Parse ATOM line and rewrite column 8 (residue name)
            # MOL2 ATOM format: index name x y z type subst_id subst_name charge
            # Columns are space-separated, need to find the subst_name field
            parts = line.split()
            if len(parts) >= 8:
                # Keep everything except replace field 7 (0-indexed) with new resname
                parts[7] = new_resname_padded.rstrip()
                output_lines.append(" ".join(parts) + "\n")
            else:
                output_lines.append(line)
        else:
            output_lines.append(line)

    with open(dst, "w") as f:
        f.writelines(output_lines)


def preprocess_pdb_for_ph(src: Path, dst: Path, ph_config: dict) -> None:
    """Prepare pH-adjusted PDB for pdb2gmx by renaming residues.

    Args:
        src: Input PDB file
        dst: Output PDB file
        ph_config: Dict with 'his_method' ('auto' or 'HIP') and 'overrides' dict
                   (keyed by "chain:resnum")
    """
    src = Path(src)
    dst = Path(dst)

    his_method = ph_config.get("his_method", "auto")
    overrides = ph_config.get("overrides", {})

    with open(src) as f:
        lines = f.readlines()

    output_lines = []

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # PDB fixed-width format:
            # columns 17-20 (0-indexed 16-19) contain the residue name
            # columns 22-26 (0-indexed 21-25) contain the residue number
            # column 22 (0-indexed 21) is the chain ID

            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = line[22:26].strip()

            # Check for per-residue override
            override_key = f"{chain}:{resnum}"
            if override_key in overrides:
                new_resname = overrides[override_key]
            elif his_method == "HIP" and resname == "HIS":
                # Global HIS->HIP rename for pH 4.9
                new_resname = "HIP"
            else:
                new_resname = resname

            # Rewrite residue name field (columns 17-20, 0-indexed 17-19)
            new_resname_padded = new_resname.rjust(3)  # Right-justify in 3-char field
            line = line[:17] + new_resname_padded + line[20:]

        output_lines.append(line)

    with open(dst, "w") as f:
        f.writelines(output_lines)


def validate_files_exist(paths: list[Path], context: str) -> None:
    """Assert all files exist.

    Args:
        paths: List of Path objects to check
        context: Description for error message

    Raises:
        FileNotFoundError: If any file is missing
    """
    missing = [p for p in paths if not p.exists()]
    if missing:
        msg = f"{context}: missing files:\n"
        for p in missing:
            msg += f"  {p}\n"
        raise FileNotFoundError(msg)


def conda_cmd(args: list[str], conda_env: str) -> list[str]:
    """Prepend conda run wrapper to command args.

    Args:
        args: Command arguments (e.g., ['antechamber', '-i', 'file.mol2', ...])
        conda_env: Conda environment name

    Returns:
        Full command as ['conda', 'run', '-n', env, '--no-capture-output', ...args]
    """
    return ["conda", "run", "-n", conda_env, "--no-capture-output"] + args
