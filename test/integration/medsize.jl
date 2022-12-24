const appfolder = "examples/CapacitatedVehicleRouting"
include("../../$appfolder/src/run.jl")

function run_mediumsize_instances_tests()
    @testset "Solving selected medium-size CVRP instance" begin
        # result = main(["../$appfolder/data/toy.vrp", "-m", "3", "-M", "3", "-u", "274.01"])
        # @test result ≈ 274.0
        result = main(["../$appfolder/data/A/A-n37-k6.vrp", "-m", "6", "-M", "6", "-u", "949.01"])
        @test result ≈ 949.0
        # result = main(["../$appfolder/data/F/F-n72-k4.vrp", "-m", "4", "-M", "4", "-u", "237.01"])
        # @test result ≈ 237.0
    end
end
