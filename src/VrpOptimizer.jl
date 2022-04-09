mutable struct VrpOptimizer
    model::VrpModel
    status_dict::Dict{MathOptInterface.TerminationStatusCode, Symbol}
end

function VrpOptimizer(model::VrpModel, _::String, _::String)
    # add the mapped arc-variables to the model
    @axis(VrpGraphs, 1:length(model.graphs))
    @variable(
        model.formulation,
        __arc[g in VrpGraphs, a in collect(keys(model.graphs[g].mappings))] >= 0
    )

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

    # add the mapping constraints to the model and inform the RCSP solver (TODO)
    # keep a Dict of VariableRefs to ids
    # use "obj = objective_function(m, AffExpr)" and "coefficien(obj, x)" to get
    # the variable costs (also used when they are mapped)
    @constraint(
        model.formulation, [v in collect(keys(map_inverse))],
        v == sum(__arc[g, a] for (g, a) in map_inverse[v])
    )

    # apply the decomposition, store the axis, create the optimizer and return it
    @dantzig_wolfe_decomposition(model.formulation, decomp, VrpGraphs)
    model.bd_graphs = decomp
    sd = Dict(
        MathOptInterface.OPTIMAL => :Optimal,
        MathOptInterface.INFEASIBLE => :Infeasible,
        MathOptInterface.ITERATION_LIMIT => :UserLimit,
        MathOptInterface.TIME_LIMIT => :UserLimit,
        MathOptInterface.OPTIMIZE_NOT_CALLED => :NotSolved
    )
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
