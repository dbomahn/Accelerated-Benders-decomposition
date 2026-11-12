# Accelerated Benders decomposition

Julia implementation accompanying the paper **“Accelerated Benders decomposition for bi-objective mixed integer programming”**.

The repository contains the supply chain network design (SCND) instance data, the data notation guide, the model and Benders implementation, and a single experiment driver.

## Repository contents

- `src/`: data reader, SCND model, Benders structures and algorithm functions.
- `scripts/DichoBD.jl`: execution script for the default accelerated Benders experiment.
- `scripts/paths.jl`: central path configuration.
- `scripts/smscndgenerator.jl`: SCND instance generator.
- `data/instances/`: supplied `.dat` instances and `Notations.txt`.
- `results/`: generated CSV output (not committed).

## Data

The supplied instance directory contains the following instances:

`test01S1.dat`, `test01S2.dat`, `test01S3.dat`, `test01S4.dat`,  
`test02S1.dat`, `test02S2.dat`, `test02S3.dat`, `test02S4.dat`,  
`test03S2.dat`, and `test03S3.dat`.

The accompanying [Notations.txt](data/instances/Notations.txt) explains the notation used for the SCND data. In broad terms, the instances describe a four-echelon supply chain consisting of suppliers, plants, distribution centres (DCs) and customers. The data cover:

- network dimensions and facility-opening limits;
- available transportation modes between consecutive echelons;
- customer demand and supplier, plant and DC capacities;
- fixed facility costs and transportation costs;
- procurement, production and warehousing costs;
- transportation, procurement, production and warehousing emissions;
- mode-specific capacity and fixed-cost data.

Each `.dat` file is read by `Data` in `src/Data.jl`. It begins with the six network-size and facility-limit values, followed by labelled data blocks. Keep the labels and block ordering intact when preparing additional instances. See [data/instances/README.md](data/instances/README.md) for further details.

## Requirements

- Julia and the packages imported by `src/AcceleratedBD.jl`.
- A working CPLEX installation and valid licence. The optimisation models use CPLEX through JuMP.

`Project.toml` defines the lightweight environment used by the automated data-reader tests. It does not yet declare or lock all optimisation and reporting dependencies used by the full experiment driver; those still need to be installed in the local Julia environment. A CPLEX installation and licence are required to solve optimisation models.

## Automated checks

GitHub Actions runs the data-reader test suite on pushes, pull requests targeting `main`, and manual workflow dispatches. The tests load every supplied `.dat` instance with `Data` from `src/Data.jl` and check parsed dimensions, key array shapes and basic data validity. They do not instantiate CPLEX, build optimisation models, or execute the Benders algorithm. Consequently, a passing CI run validates instance parsing, not solver-backed algorithmic correctness.

Run the same checks locally from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
```

## Run the default experiment

From the repository root:

```sh
julia scripts/DichoBD.jl
```

By default, the driver runs **one full accelerated Benders decomposition experiment** on `data/instances/test01S1.dat`, using:

| Parameter | Default |
|---|---:|
| Cut parameter (`pcut`) | 40 |
| Weight precision (`wdg`) | 4 |
| Weight rounding (`rort`) | `"t"` (truncate) |
| Cut coefficient precision (`dg`) | 1 |
| Buffer size (`buf`) | 10,000 |
| Objective precision (`num`) | 5 |

These defaults follow the first instance and parameter settings selected in the previous experiment driver. The script runs the non-truncated `AcceleratedBD` routine; it does not run the classical or truncated variants, repeat the experiment, or create plots.

To run a different instance, set `SCND_INSTANCE_FILE` to its path:

```sh
SCND_INSTANCE_FILE=data/instances/test02S1.dat julia scripts/DichoBD.jl
```

## Results

The driver writes four CSV files to `results/`, named using the instance filename:

- `*_pareto_front.csv`: non-dominated objective vectors;
- `*_iterations.csv`: Benders iteration log;
- `*_points.csv`: point-level algorithm log;
- `*_run_summary.csv`: experiment parameters, elapsed time and Pareto-point count.

CSV files are used consistently; they are not Excel workbooks. Existing files with the same names are overwritten by a new run.

## Path configuration

Paths can be overridden through environment variables:

- `SCND_INSTANCE_DIR`
- `SCND_RESULTS_DIR`
- `SCND_NOTATIONS_FILE`
- `SCND_GENERATOR_INPUT`
- `SCND_GENERATOR_OUTPUT_DIR`
- `SCND_GENERATOR_METADATA_FILE`
- `SCND_CLASSIC_NDP_DIR`
- `SCND_ACCELERATED_NDP_DIR`

See [docs/architecture.md](docs/architecture.md) for the path configuration used by the scripts.
