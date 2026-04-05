"""

Converts the ideal CITRIC ACID sdf file into citrate- and citrate2- species.

Citrate pkas:

Central COOH pKa = 3.13
Terminal COOH pKa = 4.76
2nd terminal COOH pKa = 6.40

As the molecule is symmetric, the two terminal COOH groups are equivalent.

"""

from pathlib import Path
import subprocess

from rdkit import Chem
from rdkit.Chem import SDWriter


REPO_ROOT = Path(__file__).parent.parent
SDF_INPUT = REPO_ROOT / "templates" / "CIT_ideal.sdf"


def find_cooh_groups(mol):
    """Return list of (carboxyl_C_idx, carbonyl_O_idx, hydroxyl_O_idx, h_idx)."""
    cooh_pattern = Chem.MolFromSmarts('[CX3:1](=[OX1:2])[OX2H1:3]')
    groups = []
    for match in mol.GetSubstructMatches(cooh_pattern):
        c_idx, o_dbl_idx, o_oh_idx = match
        o_atom = mol.GetAtomWithIdx(o_oh_idx)
        h_idx = next(
            n.GetIdx() for n in o_atom.GetNeighbors() if n.GetAtomicNum() == 1
        )
        groups.append((c_idx, o_dbl_idx, o_oh_idx, h_idx))
    return groups


def classify_cooh_groups(mol, groups):
    """
    Split COOH groups into (central, terminals).
    The central COOH carboxyl C is a direct neighbor of the quaternary C-OH carbon
    (the carbon bearing the hydroxyl and two CH2 groups).
    """
    central_c_pattern = Chem.MolFromSmarts('[CX4:1]([OX2H1])([CX4H2])[CX4H2]')
    match = mol.GetSubstructMatch(central_c_pattern)
    if not match:
        raise RuntimeError("Could not identify central quaternary carbon in citric acid")

    central_c_idx = match[0]
    central_c_neighbors = {
        n.GetIdx() for n in mol.GetAtomWithIdx(central_c_idx).GetNeighbors()
    }

    central = None
    terminals = []
    for group in groups:
        c_idx = group[0]
        if c_idx in central_c_neighbors:
            central = group
        else:
            terminals.append(group)

    if central is None:
        raise RuntimeError("Could not identify central COOH group")
    if len(terminals) != 2:
        raise RuntimeError(f"Expected 2 terminal COOH groups, found {len(terminals)}")

    return central, terminals


def apply_deprotonations(mol, deprot_pairs):
    """
    deprot_pairs: list of (h_idx, o_idx)

    Sets formal charge -1 on each specified O, then removes the H atoms
    from highest index to lowest to avoid index shifting.
    """
    rw = Chem.RWMol(mol)

    # Set all formal charges before any atom removal
    for h_idx, o_idx in deprot_pairs:
        rw.GetAtomWithIdx(o_idx).SetFormalCharge(-1)

    # Remove H atoms from highest index to lowest
    for h_idx in sorted((h for h, _ in deprot_pairs), reverse=True):
        rw.RemoveAtom(h_idx)

    Chem.SanitizeMol(rw)
    return rw.GetMol()


def sdf_to_mol2(sdf_path: Path, mol2_path: Path):
    """
    Convert SDF to MOL2 via OpenBabel.
    No --gen3d: preserves original 3D coordinates from the SDF.
    OpenBabel assigns SYBYL type O.co2 to both oxygens of COO- groups.
    """
    subprocess.run(
        ["obabel", "-i", "sdf", str(sdf_path), "-o", "mol2", "-O", str(mol2_path)],
        check=True,
    )


def main():
    mol = Chem.MolFromMolFile(str(SDF_INPUT), removeHs=False)
    if mol is None:
        raise RuntimeError(f"Failed to load {SDF_INPUT}")

    groups = find_cooh_groups(mol)
    if len(groups) != 3:
        raise RuntimeError(f"Expected 3 COOH groups in citric acid, found {len(groups)}")

    central, terminals = classify_cooh_groups(mol, groups)
    _, _, o_central_idx, h_central_idx = central
    _, _, o_term_idx, h_term_idx = terminals[0]  # terminals are equivalent by symmetry

    # citrate-: central COOH deprotonated (pKa 3.13, ~100% at pH 4.9)
    mol_minus = apply_deprotonations(mol, [(h_central_idx, o_central_idx)])

    # citrate2-: central + one terminal COOH deprotonated (terminal pKa 4.76, ~61% at pH 4.9)
    mol_2minus = apply_deprotonations(mol, [
        (h_central_idx, o_central_idx),
        (h_term_idx, o_term_idx),
    ])

    templates_dir = REPO_ROOT / "templates"
    sdf_minus = templates_dir / "CIT_minus_ideal.sdf"
    sdf_2minus = templates_dir / "CIT2_minus_ideal.sdf"

    for mol_out, path, label in [
        (mol_minus, sdf_minus, "citrate-"),
        (mol_2minus, sdf_2minus, "citrate2-"),
    ]:
        w = SDWriter(str(path))
        w.write(mol_out)
        w.close()
        print(f"Written {label:10s}: {path}  ({mol_out.GetNumAtoms()} atoms)")

    mol2_minus = REPO_ROOT / "templates" / "mol2" / "citrate_minus.mol2"
    mol2_2minus = REPO_ROOT / "templates" / "mol2" / "citrate2_minus.mol2"

    sdf_to_mol2(sdf_minus, mol2_minus)
    print(f"Written citrate-   MOL2: {mol2_minus}")
    sdf_to_mol2(sdf_2minus, mol2_2minus)
    print(f"Written citrate2-  MOL2: {mol2_2minus}")


if __name__ == "__main__":
    main()
