const appfolder = "examples/CapacitatedVehicleRouting"
include("../../$appfolder/src/run.jl")

function run_mediumsize_instances_tests()
    @testset "Solving selected medium-size CVRP instance" begin
        result = main(["../$appfolder/data/F/F-n72-k4.vrp", "-m", "4", "-M", "4", "-u", "238"])
        @test result ≈ 237.0
    end
end
