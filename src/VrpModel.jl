# Data structure to avoid cut coefficient repetitions minimizing allocations
mutable struct CutCoeffManager
    has_coeff::Vector{Int}
    nb_cuts::Int
end

CutCoeffManager() = CutCoeffManager(Int[], 0)

function nextcut!(ccm::CutCoeffManager, model::T) where {T <: AbstractVrpModel}
    ccm.nb_cuts += 1
    if isempty(model.coeffmanager.has_coeff)
        model.coeffmanager.has_coeff = zeros(Int, get_maxvarid(model))
    end
end

hascoeff(ccm::CutCoeffManager, v::Int) = (ccm.has_coeff[v] == ccm.nb_cuts)

function regcoeff!(ccm::CutCoeffManager, v::Int)
    ccm.has_coeff[v] = ccm.nb_cuts
    return
end

Coluna.@with_kw mutable struct RedCostFixAndEnumAlgorithm <:
    Coluna.Algorithm.AbstractOptimizationAlgorithm
    func::Function
end

Coluna.@with_kw mutable struct SolveByMipAlgorithm <:
    Coluna.Algorithm.AbstractOptimizationAlgorithm
    func::Function
end

mutable struct VrpModel <: AbstractVrpModel
    formulation::JuMP.Model
    form_obj::JuMP.AffExpr
    rcsp_instances::Vector{RCSPProblem}
    bd_graphs::Vector{BlockDecomposition.Root{:VrpGraphs, Int64}}
    rcc_separators::Vector{Ptr{Cvoid}}
    rank1cut_separator::Ptr{Cvoid}
    rcc_pre_separators::Vector{RCCPreSeparator}
    coeffmanager::CutCoeffManager
    variables_by_id::Vector{VariableRef}
    varids_by_var::Dict{VariableRef, Int}
    redcostfix_enum_algo::RedCostFixAndEnumAlgorithm
    solve_by_mip_algo::SolveByMipAlgorithm
    parameters::Vector{VrpParameters}
end

get_parameter(model::VrpModel, param::Symbol) =
    getfield(model.parameters[1].coluna_vrp_params, param)

get_maxvarid(model::VrpModel) = length(model.variables_by_id)

getvar(model::VrpModel, id::Int) = model.variables_by_id[id]

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
    colgenstages = Coluna.Algorithm.ColumnGeneration[]
    for stage in 1:3
        push!(colgenstages, Coluna.Algorithm.ColumnGeneration(
            pricing_prob_solve_alg = Coluna.Algorithm.SolveIpForm(optimizer_id = stage),
            smoothing_stabilization = 1.0
        ))
    end
    dummyfunc() = nothing
    redcostfix_enum_algo = RedCostFixAndEnumAlgorithm(func = dummyfunc)
    solve_by_mip_algo = SolveByMipAlgorithm(func = dummyfunc)
    colcutgen = Coluna.Algorithm.ColCutGenConquer(
        stages = colgenstages,
        primal_heuristics = [],
        before_cutgen_user_algorithm = Coluna.Algorithm.BeforeCutGenAlgo(
            redcostfix_enum_algo, "Reduced cost fixing and enumeration"
        ),
        node_finalizer = Coluna.Algorithm.NodeFinalizer(
            solve_by_mip_algo, 0, "Solver by MIP"
        )
    )
    branching = Coluna.Algorithm.StrongBranching()
    prodscore = Coluna.Algorithm.ProductScore()
    push!(branching.phases, Coluna.Algorithm.BranchingPhase(
        20, Coluna.Algorithm.RestrMasterLPConquer(), prodscore)
    )
    push!(branching.phases, Coluna.Algorithm.BranchingPhase(1, colcutgen, prodscore))
    push!(branching.rules, Coluna.Algorithm.PrioritisedBranchingRule(
        Coluna.Algorithm.SingleVarBranchingRule(), 1.0, 1.0)
    )

    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
            solver = Coluna.Algorithm.TreeSearchAlgorithm(
                branchingtreefile = "BaPTree.dot",
                conqueralg = colcutgen,
                dividealg = branching
            )
        ),
        "default_optimizer" => Gurobi.Optimizer # for the master & the subproblems
    )
    form = BlockModel(coluna) # , direct_model = true)

    # Return the VrpSolver model containing the Coluna and RCSP models
    return VrpModel(
        form, AffExpr(), RCSPProblem[],
        Vector{BlockDecomposition.Root{:VrpGraphs, Int64}}(undef, 1),
        Ptr{Cvoid}[], Ptr{Cvoid}(), RCCPreSeparator[], CutCoeffManager(), VariableRef[],
        Dict{VariableRef, Int64}(), redcostfix_enum_algo, solve_by_mip_algo, VrpParameters[]
    )
end

function set_branching_priority!(::VrpModel, ::String, ::Int)
    @warn "set_branching_priority! is not implemented... ignoring it."
end

function add_cut_callback!(::VrpModel, ::Function, ::String)
    @warn "add_cut_callback! is not implemented... ignoring it."
end
