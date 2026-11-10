#######################  small instance Benders Decomposition  ########################
using CSV,DelimitedFiles,Distributions,CPLEX,MathOptInterface,JuMP,MathOptInterface,LinearAlgebra,CPUTime,DataFrames,JLD2,vOptGeneric,StatsBase
import MathOptInterface as MOI

using PlotlyJS,DataFrames,DelimitedFiles,Colors,CSV,JLD2
#############################        2D plot      ###########################
layout = Layout(
    title="Test",
    xaxis_title="Cost",
    yaxis_title="CO2 emission",
    legend_title="Legend Title",
    font=attr(
        family="Courier New, monospace",
        size=18
    )
)


# % $CO_2$ emission unit has changed from kg to g ($\times$ 10^3) 
#######################  small instance Benders Decomposition  ########################

include(joinpath(@__DIR__, "Data.jl"))
include(joinpath(@__DIR__, "SCNDModel.jl"))
include(joinpath(@__DIR__, "BendersStructures.jl"))
include(joinpath(@__DIR__, "BendersFunctions.jl"))
