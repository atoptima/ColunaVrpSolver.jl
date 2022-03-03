mutable struct VrpSolver
    formulation::JuMP.Model
    rcsp_data::Ref{Ptr{Cvoid}}
end

function VrpSolver()
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
                branchingtreefile = "ColunaTree.dot",
                conqueralg = colgen,
                dividealg = branching
            )
        ),
        "default_optimizer" => Gurobi.Optimizer # for the master & the subproblems
    )
    form = BlockModel(coluna, direct_model = true)

    # Create an RCSP model
    rcsp = Ref{Ptr{Cvoid}}()
    ccall(
        (:createAndPrepareSolver_c, "$(ENV["RCSP_LIB_PATH"])"),
        Cvoid, (Ptr{Ptr{Cvoid}},), rcsp
    )

    # Return the VrpSolver model containing the Coluna and RCSP models
    return VrpSolver(form, rcsp)
end
