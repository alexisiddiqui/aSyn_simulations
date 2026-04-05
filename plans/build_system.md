
# System building 

When generating the system from the protein and bufffer ligands, additional ions are added to neutralize the system. 


## Flow

0. Generate topology for protein and ligands (previous step)
1. Add buffer ligands based on box size - must match stoichiometry
2. Solvate system with water
3. Add ions to condition concentration and neutralize buffer ligands


Tris+:Tris ~ 1.0:0.22 (pH 7.4)
Tris+:Tris ~ 1.0:0.11 (pH 7.2)
Citrate2-:Citrate- ~ 32:23 (pH 4.9)




# Notes

The main topology generation is separated as there will be multiple starting structures for each condition in the future.