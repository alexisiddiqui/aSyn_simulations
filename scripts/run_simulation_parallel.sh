
#!/bin/bash

set -e

MAXWORKERS=8
CONDITIONS=("tris" "saline" "extracellular" "intracellular" "lysosomal")

run_simulation() {
    local condition=$1
    echo "Starting simulation for condition: $condition"
    uv run topology_generation/run_simulation.py --n-replicates 1 --through 8_prod_1ns_test --condition "$condition" --force
    echo "Completed simulation for condition: $condition"
}

export -f run_simulation

# Run each condition in parallel with maxworkers limit
printf '%s\n' "${CONDITIONS[@]}" | xargs -P "$MAXWORKERS" -I {} bash -c 'run_simulation "$@"' _ {}



uv run topology_generation/run_simulation.py --n-replicates 1 --through 8_prod_1ns_test --condition tris --force
uv run topology_generation/run_simulation.py --n-replicates 1 --through 8_prod_1ns_test --condition saline --force
uv run topology_generation/run_simulation.py --n-replicates 1 --through 8_prod_1ns_test --condition extracellular --force
uv run topology_generation/run_simulation.py --n-replicates 1 --through 8_prod_1ns_test --condition intracellular --force
uv run topology_generation/run_simulation.py --n-replicates 1 --through 8_prod_1ns_test --condition lysosomal --force