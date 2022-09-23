const appfolder = "examples/CapacitatedVehicleRouting"
include("../../$appfolder/src/run.jl")

function run_mediumsize_instances_tests()
    @testset "Solving selected medium-size CVRP instance" begin
        main(["../$appfolder/data/F/F-n45-k4.vrp", "-m", "4", "-M", "4", "-u", "725"])
        # main(["../$appfolder/data/E/E-n51-k5.vrp", "-m", "5", "-M", "5", "-u", "522"])
        # main(["../$appfolder/data/A/A-n37-k6.vrp", "-m", "6", "-M", "6", "-u", "950"])
    end
end
