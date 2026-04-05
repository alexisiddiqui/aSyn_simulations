# Topology Generation Pipeline — aSyn Simulations

## Context

Alpha-synuclein MD simulation project needs parameterized topologies for 4 conditions (tris, extracellular, intracellular, lysosomal) before system assembly and simulation can begin. The project has MOL2 files and an aSyn PDB but lacks the GAFF2 parameterization pipeline and pdb2gmx automation. Goal: generate GROMACS-ready `.itp`/`.top` files for all ligands and pH-specific protein topologies via a config-driven pipeline with separate scripts.

---

## Files to Create

| File | Purpose |
|------|---------|
| `config/topology_config.yaml` | Shared YAML config (conditions, ligands, tool paths, protonation) |
| `topology_generation/topology_utils.py` | Shared utilities: config loading, subprocess runner, mol2/PDB fixers, validators |
| `topology_generation/generate_ligand_topology.py` | GAFF2 pipeline for tris/citrate ligands |
| `topology_generation/generate_protein_topology.py` | pdb2gmx pipeline for aSyn at each unique pH |

---

## YAML Config Schema (`config/topology_config.yaml`)

```yaml
paths:
  repo_root: "/Users/alexi/Documents/aSyn_simulations"
  gmx: "/Users/alexi/gromacs-2026/bin/gmx"
  output_dir: "topology_output"
  conda_env: "amber"

protein:
  input_pdb: "templates/protein/aSyn_s20_r1_msa1-127_n12700_do1_20260329_025853_protonated_max_plddt_12691.pdb"
  forcefield: "amber19sb"
  water_model: "opc"

conditions:
  tris:         { ph: 7.4, ligands: [TR1] }
  extracellular:{ ph: 7.4, ligands: [TR1] }
  intracellular:{ ph: 7.2, ligands: [TR1] }
  lysosomal:    { ph: 4.9, ligands: [CI1, CI2] }

ligands:
  TR1: { mol2: "templates/mol2/tris.mol2",           net_charge:  1, resname: "TR1" }
  CI1: { mol2: "templates/mol2/citrate_minus.mol2",  net_charge: -1, resname: "CI1" }
  CI2: { mol2: "templates/mol2/citrate2_minus.mol2", net_charge: -2, resname: "CI2" }

protonation:
  "7.4": { his_method: "auto", overrides: {} }
  "7.2": { his_method: "auto", overrides: {} }
  "4.9": { his_method: "HIP",  overrides: { "A:50": "HIP" } }
```

---

## Output Directory Layout

```
topology_output/
  ligands/
    TR1/  CI1/  CI2/
      {RES}_input.mol2      # resname-fixed input
      {RES}_AC.mol2         # antechamber output (GAFF2 + AM1-BCC)
      {RES}_AC.frcmod       # parmchk2 missing params
      {RES}_AC.prmtop       # tleap topology
      {RES}_AC.inpcrd       # tleap coordinates
      {RES}.amb2gmx/        # acpype working dir
      {RES}_GMX.itp         # final GROMACS params (copied here)
      {RES}_GMX.top
  protein/
    ph_7.4/  ph_7.2/  ph_4.9/
      aSyn_preprocessed.pdb
      protein.gro
      protein.top
      posre.itp
```

---

## Step 1: `topology_utils.py`

```python
def load_config(config_path: Path) -> dict: ...
def run_cmd(cmd: list[str], cwd: Path, label: str) -> None:
    # subprocess.run(check=True), write stdout/stderr to {label}.log in cwd
def fix_mol2_resname(src: Path, dst: Path, new_resname: str) -> None:
    # Plain-text rewrite of @<TRIPOS>MOLECULE name and ATOM column 8
def preprocess_pdb_for_ph(src: Path, dst: Path, ph_config: dict) -> None:
    # Rename HIS→HIP (if his_method=="HIP") and apply per-residue overrides
    # Uses PDB column 18-20 (resname) manipulation
def validate_files_exist(paths: list[Path], context: str) -> None: ...
def conda_cmd(args: list[str], conda_env: str) -> list[str]:
    # Returns ['conda', 'run', '-n', conda_env, '--no-capture-output'] + args
```

---

## Step 2: `generate_ligand_topology.py` — GAFF2 Pipeline

Collects unique ligands across all conditions, runs the full pipeline for each:

### `step1_fix_resname(mol2_src, work_dir, resname) → Path`
Calls `fix_mol2_resname()`. Required because existing mol2 files have resname `UNL1`.

### `step2_antechamber(mol2_in, work_dir, resname, net_charge, conda_env) → Path`
```
conda run -n amber --no-capture-output antechamber \
  -i {mol2_in} -fi mol2 -o {RES}_AC.mol2 -fo mol2 \
  -c bcc -nc {net_charge} -at gaff2 -rn {resname} -s 2
```

### `step3_parmchk2(ac_mol2, work_dir, resname, conda_env) → Path`
```
conda run -n amber --no-capture-output parmchk2 \
  -i {RES}_AC.mol2 -f mol2 -o {RES}_AC.frcmod -s gaff2
```

### `step4_tleap(ac_mol2, frcmod, work_dir, resname, conda_env) → (Path, Path)`
Write `tleap_{RES}.in` inline (use **absolute paths** inside tleap input):
```
verbosity 1
source leaprc.gaff2
{RES} = loadmol2 /abs/path/{RES}_AC.mol2
loadamberparams /abs/path/{RES}_AC.frcmod
check {RES}
saveamberparm {RES} /abs/path/{RES}_AC.prmtop /abs/path/{RES}_AC.inpcrd
quit
```
```
conda run -n amber --no-capture-output tleap -f tleap_{RES}.in
```

### `step5_acpype(prmtop, inpcrd, work_dir, resname) → (Path, Path)`
Uses **uv env** (not conda) — acpype is installed via uv. Use `amb2gmx` mode to skip re-running antechamber:
```
uv run acpype -p {RES}_AC.prmtop -x {RES}_AC.inpcrd -b {RES} -o gmx
```
cwd=`work_dir`. Output in `{RES}.amb2gmx/`. Copy `{RES}_GMX.itp` and `{RES}_GMX.top` to `work_dir/`.

### Validation (`validate_ligand_outputs`)
- All 6 expected files exist and non-empty
- `.itp` contains `[ moleculetype ]`, `[ atoms ]`, `[ bonds ]`
- `.frcmod` has **no** `"ATTN, NEED"` lines (missing params flag)
- Atom count in `.itp` matches source mol2
- Charge sum in `[ atoms ]` rounds to `net_charge ± 0.01 e`

---

## Step 3: `generate_protein_topology.py` — pdb2gmx Pipeline

Collects unique pH values, runs pdb2gmx once per unique pH.

### `preprocess_pdb_for_ph(src, dst, ph_config) → None`
- If `his_method == "HIP"`: rename all HIS → HIP in PDB columns 18-20
- Apply per-residue `overrides` dict (key: `"chain:resnum"`)
- PDB already has explicit H atoms — `-ignh` strips them; do not pre-remove

### `run_pdb2gmx(input_pdb, output_dir, gmx_path, forcefield, water_model, ph) → (gro, top, itp)`
```
{gmx} pdb2gmx \
  -f aSyn_preprocessed.pdb \
  -o protein.gro -p protein.top -i posre.itp \
  -ff amber19sb -water opc \
  -ignh
```
- **Non-interactive**: `-ff` and `-water` suppress all prompts
- **No `-his` flag**: HIS auto-assigned by H-bond geometry for pH 7.4/7.2; pH 4.9 handled by pre-renaming HIS→HIP
- **No `-ter` flag**: amber19sb handles N/C termini via `NMET`/`CLEU` rtp blocks automatically

### Validation (`validate_protein_outputs`)
- `protein.gro` has > 900 atom lines
- `protein.top` includes `amber19sb.ff/forcefield.itp` and `amber19sb.ff/opc.itp`
- `posre.itp` contains `[ position_restraints ]`
- pdb2gmx stderr has no lines starting with `"CRITICAL"` or `"Fatal"`

---

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| `amb2gmx` mode for acpype (not direct mol2) | Avoids acpype re-running antechamber with bundled (non-conda) AmberTools. Reads the prmtop/inpcrd from step 4. |
| Pre-rename HIS→HIP for pH 4.9 (not `-his` flag) | `-his` triggers interactive stdin prompts. Pre-renaming is simpler, deterministic, and keeps pdb2gmx fully non-interactive. |
| Inline tleap input string (not template file) | Only ~8 lines, per-ligand only by resname + abs paths. Written to disk for audit. |
| `--no-capture-output` on `conda run` | Without it, conda buffers all output until process exits — AM1-BCC can take minutes, leaving user with no feedback. |
| Separate scripts, shared config | Allows running ligand and protein steps independently. Config ensures conditions stay in sync. |
| Deduplicate by unique pH | pH 7.4 appears in tris + extracellular — generate protein topology once, reuse across conditions. |

---

## Verification

1. Run ligands: `uv run python topology_generation/generate_ligand_topology.py`
   - Check `topology_output/ligands/{TR1,CI1,CI2}/` for `_GMX.itp` files
   - Inspect `.itp` charge sums: TR1→+1, CI1→-1, CI2→-2
   - Confirm no `"ATTN, NEED"` in `.frcmod` files

2. Run protein: `uv run python topology_generation/generate_protein_topology.py`
   - Check `topology_output/protein/ph_{7.4,7.2,4.9}/` for `protein.gro`
   - Verify HIS-50 in pH 4.9 topology is `HIP` (check `.top` or `protein.gro`)
   - Validate `protein.top` includes OPC water model

3. Quick GROMACS sanity check on protein topology:
   ```
   gmx editconf -f topology_output/protein/ph_7.4/protein.gro -o test_box.gro -box 10 10 10
   ```
   Should complete without error.
