# Data structure to avoid cut coefficient repetitions minimizing allocations
mutable struct CutCoeffManager
    has_coeff::Vector{Int}
    nb_cuts::Int
end

CutCoeffManager() = CutCoeffManager(Int[], 0)

function nextcut!(ccm::CutCoeffManager, model::T) where {T <: AbstractVrpModel}
    ccm.nb_cuts += 1
    if isempty(model.coeffmanager.has_coeff)
        model.coeffmanager.has_coeff = zeros(
            Int, get_maxarcid(model.rcsp_instances[g].graph) + 1 # FIXME
        )
    end
end

hascoeff(ccm::CutCoeffManager, v::Int) = (ccm.has_coeff[v] == ccm.nb_cuts)

function regcoeff!(ccm::CutCoeffManager, v::Int)
    ccm.has_coeff[v] = ccm.nb_cuts
    return
end

mutable struct VrpModel <: AbstractVrpModel
    formulation::JuMP.Model
    rcsp_instances::Vector{RCSPProblem}
    bd_graphs::Vector{BlockDecomposition.Root{:VrpGraphs, Int64}}
    rcc_separators::Vector{Ptr{Cvoid}}
    coeffmanager::CutCoeffManager
    variables_by_id::Vector{VariableRef}
    varids_by_var::Dict{VariableRef, Int}
end

get_maxvarid(model::VrpModel) = length(model.variables)

getvar(model::VrpModel, id::Int) = model.variables[id]

function getvarid!(model::VrpModel, var::VariableRef)
    new_varid = get_maxvarid(model) + 1
    varid = get(model.varids_by_var, var, new_varid)
    if varid == new_varid
        resize!(model.variables_by_id, new_varid)
        model.variables_by_id[new_varid] = var
        model.varids_by_var[var] = new_varid
    end
    return varid
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
        Ptr{Cvoid}[], CutCoeffManager(), VariableRef[]
    )
end

function set_branching_priority!(::VrpModel, ::String, ::Int)
    @warn "set_branching_priority! is not implemented... ignoring it."
end