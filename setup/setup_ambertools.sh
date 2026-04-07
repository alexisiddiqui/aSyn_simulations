mamba create -n amber -c conda-forge ambertools -y
mamba create -n psiresp -c conda-forge psiresp -y
conda run -n psiresp pip install "pydantic<2.0"