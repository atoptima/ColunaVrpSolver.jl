mutable struct VrpOptimizer
    model::VrpModel
    status_dict::Dict{MathOptInterface.TerminationStatusCode, Symbol}
end

function VrpOptimizer(model::VrpModel, _::String, _::String)
    # add the mapped arc-variables to the model
    A = [collect(keys(model.graphs[g].mappings)) for g in 1:length(model.graphs)]
    @axis(VrpGraphs, 1:length(model.graphs))
    @variable(model.formulation, __arc[g in VrpGraphs, a in A[g]] >= 0)

    # compute the inverse of the mapping function
    map_inverse = Dict{VariableRef, Vector{Tuple{Int, Int}}}()
    for graph in model.graphs
        for (arcid, mapped_modelvars) in graph.mappings
            for modelvar in mapped_modelvars
                mapped_arcs = get(map_inverse, modelvar, Tuple{Int, Int}[])
                if isempty(mapped_arcs)
                    map_inverse[modelvar] = mapped_arcs
                end
                push!(mapped_arcs, (graph.id, arcid))
            end
        end
    end

    # add the mapping constraints to the model
    # TODO when switching to implicit master:
    # - keep a Dict of VariableRefs to ids
    # - use "obj = objective_function(m, AffExpr)" and "coefficient(obj, x)" to get
    #   the variable costs (also used when they are mapped)
    # - inform the RCSP solver about mappings and variable costs (for enumeration)
    @constraint(
        model.formulation, [v in collect(keys(map_inverse))],
        v == sum(__arc[g, a] for (g, a) in map_inverse[v])
    )

    # apply the decomposition and store the axis
    @dantzig_wolfe_decomposition(model.formulation, decomp, VrpGraphs)
    model.bd_graphs = decomp

    # preallocate a vector to store the reduced costs
    arc_rcosts = [zeros(Float64, graph.max_arcid + 1) for graph in model.graphs]

    # Define the function to perform pricing via RCSP solver
    function solve_RCSP_pricing(cbdata)

        # Get the reduced costs of arc variables
        g = BD.callback_spid(cbdata, model.formulation)
        for a in A[g]
            arc_rcosts[a + 1] = BD.callback_reduced_cost(cbdata, __arc[g, a])
        end
        @show arc_rcosts

        # TODO: ccall :runPricing_c and :getOutputPaths_c
        #= submit the routes with the smallest reduced costs
        routes_idx = [r for r in 1:length(spm.route_vars[rt])]
        nb_cols = min(length(routes_idx), max_cols)
        cols_idx = partialsort!(routes_idx, 1:nb_cols, by = (r -> rc_routes[r]))
        # @show rc_routes[routes_idx[1]]
        # println("rc_routes = $([rc_routes[r] for r in cols_idx])")
        for r in cols_idx
            solvals = Float64[]
            solvars = JuMP.VariableRef[]
            solvisits = UInt512(0)
            i = 0
            for j in spm.route_vars[rt][r].visits
                e = (i < j) ? (i, j) : (j, i)
                push!(solvals, 1.0)
                push!(solvars, x[rt, e])
                solvisits |= (UInt512(1) << (j - 1))
                i = j
            end
            push!(solvals, 1.0)
            push!(solvars, x[rt, (0, i)])
            push!(solvals, 1.0)
            push!(solvars, y[rt])
            push!(solvals, spm.route_vars[rt][r].cost)
            push!(solvars, _cost[rt])
            # @show tuple.(solvars, solvals)
            MOI.submit(
                sp, BD.PricingSolution(cbdata), rc_routes[r], solvars, solvals,
                RouteVarData(solvisits, spm.route_vars[rt][r].visits)
            )
        end =#
    end

    # create a dictionary of return values for compatibility with old applications
    sd = Dict(
        MathOptInterface.OPTIMAL => :Optimal,
        MathOptInterface.INFEASIBLE => :Infeasible,
        MathOptInterface.ITERATION_LIMIT => :UserLimit,
        MathOptInterface.TIME_LIMIT => :UserLimit,
        MathOptInterface.OPTIMIZE_NOT_CALLED => :NotSolved
    )

    # create the optimizer and return it
    return VrpOptimizer(model, sd)
end

function set_cutoff!(opt::VrpOptimizer, cutoffvalue::Float64)
    objectiveprimalbound!(opt.model.formulation, cutoffvalue)
end

function optimize!(opt::VrpOptimizer)
    optimize!(opt.model.formulation)
    status = get(opt.status_dict, JuMP.termination_status(opt.model.formulation), :Error)
    hassol = (JuMP.result_count(opt.model.formulation) > 0)
    return status, hassol
end
