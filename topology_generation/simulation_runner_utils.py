"""Utilities for running GROMACS simulation stages.

Provides:
  - run_grompp: wrapper for gmx grompp
  - run_mdrun: wrapper for gmx mdrun
  - run_make_ndx: wrapper for gmx make_ndx
  - StageConfig: dataclass for stage metadata
"""

import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class StageConfig:
    """Metadata for a single simulation stage."""

    name: str
    mdp: str
    use_posres_ref: bool
    continuation: bool
    replicated: bool
    depends_on: str = None


def run_gmx(gmx_path: str, args: list, cwd: Path, label: str, stdin: str = None) -> None:
    """Run a gmx command, capturing output to {label}.log.

    Args:
        gmx_path: Path to gmx executable
        args: Command arguments (e.g. ['grompp', '-f', 'file.mdp', ...])
        cwd: Working directory
        label: Prefix for log file
        stdin: String to pipe to stdin (e.g. "q\n" for make_ndx)

    Raises:
        RuntimeError: If gmx exits non-zero
    """
    log_file = cwd / f"{label}.log"

    cmd = [gmx_path] + args

    try:
        with open(log_file, "w") as log:
            result = subprocess.run(
                cmd,
                cwd=cwd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                input=stdin,
            )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"{label} timed out after 1 hour. Log: {log_file}")
    except Exception as e:
        raise RuntimeError(f"{label} failed: {e}. Log: {log_file}")

    if result.returncode != 0:
        raise RuntimeError(f"{label} exited with code {result.returncode}. Log: {log_file}")

    print(f"  → {label} complete ({log_file})")


def run_make_ndx(gmx_path: str, gro_file: Path, ndx_file: Path, cwd: Path, force: bool = False) -> None:
    """Run gmx make_ndx to generate index file.

    Args:
        gmx_path: Path to gmx executable
        gro_file: Input coordinate file
        ndx_file: Output index file
        cwd: Working directory
        force: If True, overwrite exist file
    """
    if ndx_file.exists() and not force:
        print(f"  Index file already exists: {ndx_file.name}")
        return

    args = ["make_ndx", "-f", str(gro_file), "-o", str(ndx_file)]
    run_gmx(gmx_path, args, cwd, "make_ndx", stdin="q\n")


def run_grompp(
    gmx_path: str,
    mdp_file: Path,
    coord_file: Path,
    top_file: Path,
    ndx_file: Path,
    tpr_file: Path,
    cwd: Path,
    label: str,
    ref_file: Path = None,
    cpt_file: Path = None,
) -> None:
    """Run gmx grompp.

    Args:
        gmx_path: Path to gmx executable
        mdp_file: MDP input file
        coord_file: Coordinate input file (-c)
        top_file: Topology file (-p)
        ndx_file: Index file (-n)
        tpr_file: TPR output file (-o)
        cwd: Working directory
        label: Log label
        ref_file: Reference coords for restraints (-r), optional
        cpt_file: Checkpoint for continuation (-t), optional
    """
    args = [
        "grompp",
        "-f",
        str(mdp_file),
        "-c",
        str(coord_file),
        "-p",
        str(top_file),
        "-n",
        str(ndx_file),
        "-o",
        str(tpr_file),
        "-maxwarn",
        "2",
    ]

    if ref_file and ref_file.exists():
        args.extend(["-r", str(ref_file)])

    if cpt_file and cpt_file.exists():
        args.extend(["-t", str(cpt_file)])

    run_gmx(gmx_path, args, cwd, label)


def run_mdrun(gmx_path: str, tpr_file: Path, cwd: Path, label: str = "mdrun") -> None:
    """Run gmx mdrun.

    Args:
        gmx_path: Path to gmx executable
        tpr_file: TPR input file (-s)
        cwd: Working directory
        label: Log label
    """
    args = [
        "mdrun",
        "-v",
        "-deffnm",
        "md",
        "-s",
        str(tpr_file),
        "-ntmpi",
        "1",
        "-ntomp",
        "8",
        "-pin",
        "on",
    ]
    run_gmx(gmx_path, args, cwd, label)
