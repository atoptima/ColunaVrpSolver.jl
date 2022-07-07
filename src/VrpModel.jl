# Data structure to avoid cut coefficient repetitions minimizing allocations
mutable struct CutCoeffManager
    has_coeff::Vector{Vector{Int}}
    nb_cuts::Int
end

CutCoeffManager() = CutCoeffManager(Vector{Int}[], 0)

function nextcut!(ccm::CutCoeffManager, model::T) where {T <: AbstractVrpModel}
    ccm.nb_cuts += 1
    if isempty(model.coeffmanager.has_coeff)
        resize!(model.coeffmanager.has_coeff, length(model.rcsp_instances))
        for g in 1:length(model.rcsp_instances)
            model.coeffmanager.has_coeff[g] = zeros(
                Int, model.rcsp_instances[g].graph.max_arcid + 1
            )
        end
    end
end

hascoeff(ccm::CutCoeffManager, g::Int, a::Int) = (ccm.has_coeff[g][a + 1] == ccm.nb_cuts)

function regcoeff!(ccm::CutCoeffManager, g::Int, a::Int)
    ccm.has_coeff[g][a + 1] = ccm.nb_cuts
    return
end

mutable struct VrpModel <: AbstractVrpModel
    formulation::JuMP.Model
    rcsp_instances::Vector{RCSPProblem}
    bd_graphs::Vector{BlockDecomposition.Root{:VrpGraphs, Int64}}
    rcc_separators::Vector{Ptr{Cvoid}}
    coeffmanager::CutCoeffManager
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
    return VrpModel(
        form, RCSPProblem[],
        Vector{BlockDecomposition.Root{:VrpGraphs, Int64}}(undef, 1),
        Ptr{Cvoid}[], CutCoeffManager()
    )
end

function set_branching_priority!(::VrpModel, ::String, ::Int)
    @warn "set_branching_priority! is not implemented... ignoring it."
end