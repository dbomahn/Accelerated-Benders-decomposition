# Run one reproducible accelerated Benders experiment.
#
# Defaults reproduce the first experiment selected by the original driver:
# test01S1.dat, pcut = 40, wdg = 4, rort = "t", dg = 1, buf = 10_000.
# Override SCND_INSTANCE_FILE to use another instance file.

include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(@__DIR__, "..", "src", "AcceleratedBD.jl"))

using CSV
using DataFrames

# ---------------------------------------------------------------------------
# Experiment configuration
# ---------------------------------------------------------------------------

const INSTANCE_FILE = get(
    ENV,
    "SCND_INSTANCE_FILE",
    joinpath(INSTANCE_DIR, "test01S1.dat"),
)

const PCUT = 40
const WDG = 4
const RORT = "t"
const DG = 1
const BUFFER_SIZE = 10_000

# The Benders implementation reads this precision as a global variable.
global num = 5

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

isfile(INSTANCE_FILE) || error(
    "Instance file not found: $(INSTANCE_FILE). " *
    "Set SCND_INSTANCE_FILE to the path of a .dat instance file.",
)

mkpath(RESULTS_DIR)

instance = Data(INSTANCE_FILE)
instance_name = splitext(basename(INSTANCE_FILE))[1]

println("Running accelerated Benders decomposition")
println("Instance: $(INSTANCE_FILE)")
println("Parameters: pcut=$(PCUT), wdg=$(WDG), rort=$(RORT), dg=$(DG), buffer=$(BUFFER_SIZE)")

elapsed = @CPUelapsed begin
    solution_points, iteration_log, point_log = AcceleratedBD(
        instance,
        PCUT,
        WDG,
        RORT,
        DG,
        BUFFER_SIZE,
    )
end

pareto_front = NDfilter(solution_points)
pareto_results = DataFrame(
    objective_1 = [point[1] for point in pareto_front],
    objective_2 = [point[2] for point in pareto_front],
)

CSV.write(
    joinpath(RESULTS_DIR, "$(instance_name)_pareto_front.csv"),
    pareto_results,
)
CSV.write(
    joinpath(RESULTS_DIR, "$(instance_name)_iterations.csv"),
    iteration_log,
)
CSV.write(
    joinpath(RESULTS_DIR, "$(instance_name)_points.csv"),
    point_log,
)

run_summary = DataFrame(
    instance = [basename(INSTANCE_FILE)],
    elapsed_seconds = [round(elapsed; digits = 2)],
    pcut = [PCUT],
    wdg = [WDG],
    rort = [RORT],
    dg = [DG],
    buffer_size = [BUFFER_SIZE],
    precision_digits = [num],
    pareto_points = [nrow(pareto_results)],
)
CSV.write(
    joinpath(RESULTS_DIR, "$(instance_name)_run_summary.csv"),
    run_summary,
)

println("Completed in $(round(elapsed; digits = 2)) seconds.")
println("Pareto points: $(nrow(pareto_results))")
println("Results written to: $(RESULTS_DIR)")
