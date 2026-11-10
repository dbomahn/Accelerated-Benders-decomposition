# SCND instance data

This directory contains the supplied SCND instance files (`.dat`) and the notation guide, [Notations.txt](Notations.txt).

## Instance files

- `test01S1.dat` to `test01S4.dat`
- `test02S1.dat` to `test02S4.dat`
- `test03S2.dat` and `test03S3.dat`

The files encode four-echelon supply chain network data: suppliers, plants, distribution centres (DCs) and customers. The reader in `src/Data.jl` expects the six leading values (network dimensions and facility limits), followed by labelled blocks in the order used by the supplied files. Do not remove or rename those block labels.

## Notation guide

`Notations.txt` describes the network dimensions, transportation modes, demand and capacities, fixed and variable costs, emissions, mode capacities, mode fixed costs and additional network attributes. Some entries describe broader generator or network data; the optimisation data reader consumes the fields explicitly parsed in `src/Data.jl`.

## Selecting an instance

The default experiment uses `test01S1.dat`. To select another instance, set `SCND_INSTANCE_FILE` to its path when running `scripts/DichoBD.jl`. To relocate the instance directory, set `SCND_INSTANCE_DIR` as described in the main README.
