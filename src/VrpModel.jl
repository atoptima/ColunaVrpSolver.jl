mutable struct VrpModel
    formulation::JuMP.Model
    graphs::Vector{VrpGraph}
    rcsp_instances::Vector{Ptr{Cvoid}}
end

function VrpModel()
    # Create a Coluna model
    colgen = Coluna.Algorithm.ColCutGenConquer(
        stages = [Coluna.Algorithm.ColumnGeneration()],
        primal_heuristics = [],
    )
    branching = Coluna.Algorithm.StrongBranching()
    push!(branching.phases, Coluna.Algorithm.BranchingPhase(
        20, Coluna.Algorithm.RestrMasterLPConquer())
    )
    push!(branching.phases, Coluna.Algorithm.BranchingPhase(1, colgen))
    push!(branching.rules, Coluna.Algorithm.PrioritisedBranchingRule(
        Coluna.Algorithm.SingleVarBranchingRule(), 1.0, 1.0)
    )

    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
            solver = Coluna.Algorithm.TreeSearchAlgorithm(
                branchingtreefile = "BaPTree.dot",
                conqueralg = colgen,
                dividealg = branching
            )
        ),
        "default_optimizer" => Gurobi.Optimizer # for the master & the subproblems
    )
    form = BlockModel(coluna, direct_model = true)

    # Return the VrpSolver model containing the Coluna and RCSP models
    return VrpModel(form, BlockDecomposition.Axis(:noname, []), Ptr{Cvoid}[])
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

    # apply the decomposition and return it
    @dantzig_wolfe_decomposition(model.formulation, decomp, VrpGraphs)
    return decomp
end