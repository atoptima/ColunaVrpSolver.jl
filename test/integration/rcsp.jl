@testset "Calling the C++ library to solve the RCSP problem" begin
    @test haskey(ENV, "RCSP_LIB_PATH")

    # build a toy VRP model with one graph
    toy = VrpModel()
    V = [0, 1, 2, 3, 4]
    d(i) = (i == 0) ? 0 : (i + 1)
    G = VrpGraph(toy, V, 0, 0, (0, 4))
    resid = add_resource!(G, main = true)
    for i in V
        set_resource_bounds!(G, i, resid, 0.0, 10.0)
    end
    max_arcid = 0
    id_to_arc = Dict{Int,Tuple{Int, Int}}()
    for i in 0:3
        for j in (i + 1):4
            arcid = add_arc!(G, i, j)
            id_to_arc[arcid] = (i, j)
            max_arcid = max(max_arcid, arcid)
            set_arc_consumption!(G, arcid, resid, (d(i) + d(j)) / 2)
            arcid = add_arc!(G, j, i)
            id_to_arc[arcid] = (j, i)
            max_arcid = max(max_arcid, arcid)
            set_arc_consumption!(G, arcid, resid, (d(i) + d(j)) / 2)
        end
    end
    add_graph!(toy, G)
    rcsp = toy.rcsp_instances[1] 
    @test rcsp.solver != Ptr{Cvoid}(0)

    # set the reduced costs of variables assuming one-to-one correspondence
    # between arc ids and variable ids (valid only for explicit master)
    var_rcosts = zeros(Float64, max_arcid + 1)
    for (arcid, (i, j)) in id_to_arc
        var_rcosts[arcid + 1] = abs(j - i) - 100.0 * d(j)
    end
    paths = ColunaVrpSolver.run_rcsp_pricing(rcsp, 0, var_rcosts)
    @test paths == [
        [0, 8, 9, 8, 3], [2, 9, 8, 9, 1], [2, 14, 15, 3], [6, 17, 9, 1],
        [2, 16, 13, 1], [0, 12, 17, 3], [0, 8, 16, 7], [6, 13, 8, 3],
        [4, 11, 10, 5], [2, 9, 12, 7], [4, 15, 9, 1], [0, 10, 15, 3],
        [2, 14, 11, 1], [0, 8, 14, 5], [4, 11, 8, 3], [6, 19, 5],
        [2, 9, 10, 5], [4, 18, 7], [0, 12, 13, 1], [0, 10, 11, 1],
        [2, 9, 8, 3], [6, 17, 3], [2, 16, 7], [0, 8, 9, 1],
        [4, 15, 3], [2, 14, 5], [6, 13, 1], [0, 12, 7],
        [0, 10, 5], [4, 11, 1], [2, 9, 1], [0, 8, 3],
        [6, 7], [4, 5], [2, 3], [0, 1]
    ]
end
