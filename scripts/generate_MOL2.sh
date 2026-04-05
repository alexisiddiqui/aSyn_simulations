# first generate tris.mol2 using openbabel from templates/TRS_ideal.sdf


obabel -i sdf templates/TRS_ideal.sdf -o mol2 -O templates/mol2/tris.mol2 --gen3d
