mutable struct VrpModel <: AbstractVrpModel
    formulation::JuMP.Model
    rcsp_instances::Vector{RCSPProblem}
    bd_graphs::BlockDecomposition.Axis
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
    return VrpModel(form, RCSPProblem[], BlockDecomposition.Axis(:noname, []))
end
