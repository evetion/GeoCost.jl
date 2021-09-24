using GeoCost
using Test
using BenchmarkTools

include("pcraster.jl")


const nodata = NaN

@testset "GeoCost.jl" begin
    @testset "spread" begin

        # Examples on https://pcraster.geo.uu.nl/pcraster/4.3.0/documentation/pcraster_manual/sphinx/op_spread.html#spread

        result = [6.83 205 6.83 6 4; 204 4 4.83 9.83 9; 2.83 2 7 11.9 17; 2 0 4 207 221; 2.83 NaN 5.66 209 413]
        points = [0.0 0 0 0 2; 0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0]
        initial = [8.0 8 8 8 4; 8 8 8 8 8; 8 8 8 8 8; 0 0 8 8 8; 0 0 8 8 8]
        friction = [1.0 200 1 1 1; 200 1 1 4 4; 1 1 4 4 4; 1 1 3 200 200; 1 NaN 3 200 4]
        # @test result == spread(points, initial, friction)

        result = [2.8284 2 2 2 0; 2 0 0 0 NaN; 2 0 2 2 2.8284; 2 0 2 4 4.8284; 2 0 2 4 6]
        points = [0.0 0 0 0 6; 0 1 1 2 NaN; 0 4 0 0 0; 0 2 0 0 0; 0 3 0 0 0]

        # PCRaster
        presult = spread_pcr(points, 0, 1)
        @test all(result[.~isnan.(result)] .≈ presult[.~isnan.(presult)])

        # Julian
        @test result == spread(points, 0.0, 1.0; res=2.0)
    end
end
