@testset "Calling the C++ library to solve the RCSP problem" begin
    @test haskey(ENV, "RCSP_LIB_PATH")
    toy = VrpModel()
    V = [0, 1, 2, 3, 4]
    d(i) = (i == 0) ? 0 : (i + 1)
    G = VrpGraph(toy, V, 0, 0, (0, 4))
    resid = add_resource!(G, main = true)
    for i in V
        set_resource_bounds!(G, i, resid, 0.0, 10.0)
    end
    for i in 0:3
        for j in (i + 1):4
            arcid = add_arc!(G, i, j)
            set_arc_consumption!(G, arcid, resid, (d(i) + d(j)) / 2)
            arcid = add_arc!(G, j, i)
            set_arc_consumption!(G, arcid, resid, (d(i) + d(j)) / 2)
        end
    end
    add_graph!(toy, G)
    @test toy.rcsp_instances[1] != Ptr{Nothing}(0)
end
