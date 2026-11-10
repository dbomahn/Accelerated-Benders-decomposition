# Accelerated Benders decomposition

Julia implementation for the paper **Accelerated Benders decomposition for
bi-objective mixed integer programming**. 

Set paths through environment variables rather than editing source files:
```sh
export SCND_INSTANCE_DIR=/path/to/instances
export SCND_NOTATIONS_FILE=/path/to/Notations.txt
export SCND_RESULTS_DIR=/path/to/results
julia scripts/DichoBD.jl
```
See [docs/architecture.md](docs/architecture.md) for all path variables.
