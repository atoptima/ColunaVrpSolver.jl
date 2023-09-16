const appfolder = "examples/CapacitatedVehicleRouting"
include("../../$appfolder/src/run.jl")

function run_mediumsize_instances_tests()
    Profile.init(; n = 10^6, delay = 0.005)

    @testset "Solving selected medium-size CVRP instance" begin
        # result = main(["../$appfolder/data/toy.vrp", "-m", "3", "-M", "3", "-u", "274.01"])
        # @test result ≈ 274.0

        # result = main([
        #     "../$appfolder/data/M/M-n151-k12.vrp",
        #     "-m", "12",
        #     "-M", "12",
        #     "-u", "1015.01",
        #     "-c", "../$appfolder/config/CVRP.cfg"
        # ])
        result = main([
            "../$appfolder/data/A/A-n37-k6.vrp",
            "-m", "6",
            "-M", "6",
            "-u", "949.01",
            "-c", "../$appfolder/config/CVRP_0.cfg"
        ])
        @test result ≈ 949.0

        # @profile main(["../$appfolder/data/A/A-n37-k6.vrp", "-m", "6", "-M", "6", "-u", "949.01"])
        # pprof(; webport = 58599)
        # println("View profile at http://localhost:58599/")
        # while true
        #     sleep(1.0)
        # end
        # result = main(["../$appfolder/data/M/M-n151-k12.vrp", "-m", "12", "-M", "12", "-u", "1015.01"])
        # @test result ≈ 1015.0
        # result = main(["../$appfolder/data/F/F-n72-k4.vrp", "-m", "4", "-M", "4", "-u", "237.01"])
        # @test result ≈ 237.0
    end
end
