# Path configuration

All local paths are configured in scripts/paths.jl and can be overridden
with environment variables: SCND_INSTANCE_DIR, SCND_RESULTS_DIR,
SCND_NOTATIONS_FILE, SCND_GENERATOR_INPUT, SCND_GENERATOR_OUTPUT_DIR,
SCND_GENERATOR_METADATA_FILE, SCND_CLASSIC_NDP_DIR and
SCND_ACCELERATED_NDP_DIR.

Example:
```sh
export SCND_INSTANCE_DIR=/path/to/instances
export SCND_NOTATIONS_FILE=/path/to/Notations.txt
export SCND_RESULTS_DIR=/path/to/results
julia scripts/DichoBD.jl
```

The generator accepts the generator-specific variables above. The original
algorithm and experiment sequence are otherwise retained. CPLEX requires a
local installation and licence.
