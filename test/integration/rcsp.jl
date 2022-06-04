using JuMP

function run_rcsp_integration_tests()

    # Define the toy CVRP instance data 
    V = [0, 1, 2, 3, 4]
    V⁺ = [1, 2, 3, 4]
    E = [(i, j) for i in 0:3 for j in (i + 1):4]
    c(e) = abs(e[1] - e[2])
    δ(i) = [(_i, _j) for (_i, _j) in E if i in (_i, _j)]
    d(i) = (i == 0) ? 0.0 : Float64(i + 1)
    Q = 10.0

    function build_toy_model(
        ; with_ngpaths::Bool = false, with_capacitycuts::Bool = false
    )
        # create the master problem
        toy = VrpModel()
        @variable(toy.formulation, x[e in E], Int)
        @objective(toy.formulation, Min, sum(c(e) * x[e] for e in E))
        @constraint(toy.formulation, deg[i in V⁺], sum(x[e] for e in δ(i)) == 2.0)
     
        # create the subproblem graph
        G = VrpGraph(toy, V, 0, 0, (0, 4))
        resid = add_resource!(G, main = true)
        for i in V
            set_resource_bounds!(G, i, resid, 0.0, Q)
        end
        max_arcid = 0
        id_to_arc = Dict{Int,Tuple{Int, Int}}()
        for (i, j) in E
            arcid = add_arc!(G, i, j)
            add_arc_var_mapping!(G, arcid, x[(i,j)])
            id_to_arc[arcid] = (i, j)
            max_arcid = max(max_arcid, arcid)
            set_arc_consumption!(G, arcid, resid, (d(i) + d(j)) / 2)
            arcid = add_arc!(G, j, i)
            add_arc_var_mapping!(G, arcid, x[(i,j)])
            id_to_arc[arcid] = (j, i)
            max_arcid = max(max_arcid, arcid)
            set_arc_consumption!(G, arcid, resid, (d(i) + d(j)) / 2)
        end
        add_graph!(toy, G)
        if with_ngpaths
            set_vertex_packing_sets!(toy, [[(G,i)] for i in V⁺])
            define_elementarity_sets_distance_matrix!(
                toy, G, [[Float64(c((i, j))) for j in V⁺] for i in V⁺]
            )
        end
        if with_capacitycuts
            add_capacity_cut_separator!(toy, [ ( [(G,i)], d(i) ) for i in V⁺], Q)
        end

        return toy, x, id_to_arc, max_arcid
    end

    @testset "Calling the C++ library to solve the RCSP problem" begin
        @test haskey(ENV, "RCSP_LIB_PATH")

        # build a toy VRP model with one graph
        toy, _, id_to_arc, max_arcid = build_toy_model()
        ColunaVrpSolver.build_solvers!(toy)
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

    @testset "Solving a complete CVRP toy instance" begin
        for (ng, cc) in [(false, false), (true, false), (true, true)]
            # build a toy CVRP model and solve
            toy, x, _, _ = build_toy_model(
                with_ngpaths = ng, with_capacitycuts = cc
            )
            if cc
                @test toy.rcc_separators[1] != Ptr{Cvoid}(0)
            end
            opt = VrpOptimizer(toy, "", "toy")
            set_cutoff!(opt, 100.0)
            (status, solution_found) = optimize!(opt)

            # print and test the result
            @show status, solution_found
            @test (status, solution_found) == (:Optimal, true)
            if solution_found
                obj = get_objective_value(opt)
                @show obj
                @test obj ≈ 12.0
                sol = Float64[]
                for e in E
                    val = get_value(opt, x[e])
                    push!(sol, val)
                    if val > 1e-5
                        println("x[$e] = $val")
                    end
                end
                @test reduce(&,
                    sol .≈ [1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
                )
            end
        end
    end
end
