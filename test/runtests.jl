using Test
using DelimitedFiles

include(joinpath(@__DIR__, "..", "src", "Data.jl"))

const INSTANCE_DIR = joinpath(@__DIR__, "..", "data", "instances")
const INSTANCE_FILES = sort(filter(
    path -> endswith(path, ".dat"),
    readdir(INSTANCE_DIR; join = true),
))

@testset "SCND instance data reader" begin
    @test !isempty(INSTANCE_FILES)

    for instance_file in INSTANCE_FILES
        @testset "$(basename(instance_file))" begin
            data = Data(instance_file)

            @test data.I > 0
            @test data.J > 0
            @test data.K > 0
            @test data.L > 0

            @test size(data.Mij) == (data.I, data.J)
            @test size(data.Mjk) == (data.J, data.K)
            @test size(data.Mkl) == (data.K, data.L)

            @test length(data.demand) == data.L
            @test length(data.cas) == data.I
            @test length(data.cap) == data.J
            @test length(data.cad) == data.K

            @test size(data.fc) == (data.J + data.K, 2)
            @test length(data.vcs) == data.I
            @test size(data.vcp) == (data.J, 2)
            @test size(data.vcd) == (data.K, 2)

            @test all(isfinite, data.demand)
            @test all(>(0), data.demand)
            @test all(isfinite, data.fc)
        end
    end
end
