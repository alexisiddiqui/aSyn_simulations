# aSyn_simulations


Simulations for aSyn in different environments.

4 conditions:
1. extracellular
2. intracellular
3. lysosomal
4. tris buffer

2 buffers: 
1. tris
2. citrate

4 different ions: 
1. Na+
2. K+
3. Ca2+
4. Mg2+




# Conditions

| ion | extracellular | intracellular | lysosomal | tris | saline
| :--- | :--- | :--- | :--- | :--- | :--- |
| Na+ | 143 | 15 | 20 | 0 | 150
| K+ | 4 | 140 | 60 | 0 | 150
| Ca2+ | 2.5 | 100 nM | 0 | 0 | 0
| Mg2+ | 0.7 | 10 | 0 | 0 | 0
| pH | 7.4 | 7.2 | 4.9 | 7.4 | 7.4
| buffer | 20 mM Tris | 20 mM Tris | 20 mM citrate | 20 mM Tris | --

- use pyyaml to handle configration for simulation setup?

# Small molecules (buffers)
- Tris+
- Citrate- and Citrate2-

Tris+:Tris ~ 1.0:0.22 (pH 7.4)
Tris+:Tris ~ 1.0:0.11 (pH 7.2)
Citrate2-:Citrate- ~ 32:23 (pH 4.9)

At a box size of 10nm we need 12 buffer molecules to keep the concentration at 20mM.

# Simulation setup
## Force fields and software
- GROMACS 26 - Amberff19 + OPC water
- Generate ligands (tris (TR0), tris+ (TR1),citrate- (CI1),citrate2- (CI2))


Small molecule (MOL2)
        ↓
  Charge assignment (AM1-BCC or RESP)
        ↓
  GAFF2 atom typing (antechamber)
        ↓
  Parameter generation (parmchk2)
        ↓
  AMBER topology (tleap)
        ↓
  Convert to GROMACS (acpype)

# Order of operations:

Tris, saline, extracellular, intracellular, lysosomal

1. First test that Tris and saline work as expected. 
    - Generate topology 
    - Short 100 ns of simulations.

2. Then proceed to Extra/Intra and then lysosomal.

3. Start with 1 starting structure (3 replicates) for each condition. Move to multiple starting structures once the protocol is validated.


Amber tools available at 'amber' conda environment.
# Setup ambertools env
conda create -n amber -c conda-forge ambertools
