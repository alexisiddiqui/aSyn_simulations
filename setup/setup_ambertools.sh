# ambertools
mamba create -n amber -c conda-forge ambertools -y

# psiresp (for RESP charge calculations on macOS Arm64)
mamba create -n psiresp python=3.10 psi4=1.10 -c conda-forge -y
conda run -n psiresp pip install psiresp "pydantic<2.0"
mamba install -n psiresp -c conda-forge rdkit -y