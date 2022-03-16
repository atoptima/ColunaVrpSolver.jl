@testset "Calling the C++ library to solve the RCSP problem" begin
    @test haskey(ENV, "RCSP_LIB_PATH")
    toy = VrpModel()
    G = VrpGraph(toy, [0, 1, 2, 3, 4], 0, 0, (0, 4))
    resid = add_resource!(G, main = true)
    add_graph!(toy, G)
    @test toy.rcsp_instances[1] != Ptr{Nothing}(0)
end
