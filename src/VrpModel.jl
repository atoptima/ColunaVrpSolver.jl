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
    cutsep_phase::Int
    rcc_pre_separators::Vector{RCCPreSeparator}
    coeffmanager::CutCoeffManager
    variables_by_id::Vector{VariableRef}
    varids_by_var::Dict{VariableRef, Int}
    spids_by_var::Dict{VariableRef, Vector{Int}}
    redcostfix_enum_algo::RedCostFixAndEnumAlgorithm
    solve_by_mip_algo::SolveByMipAlgorithm
    parameters::Vector{VrpParameters}
    branch_priors::Dict{String, Int}
    nb_subproblems::Int
    cfg_fname::String
    cutoffvalue::Float64
    packing_sets::Vector{Vector{Tuple{Int, Int}}}
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
    dummyfunc() = nothing
    redcostfix_enum_algo = RedCostFixAndEnumAlgorithm(func = dummyfunc)
    solve_by_mip_algo = SolveByMipAlgorithm(func = dummyfunc)

    # Create a Coluna model
    if rcsp_path == ""
        tree_search = BapcodTreeSearchWrapper(Coluna.Optimizer[], VrpModel[])
    else
        colgen = Coluna.Algorithm.ColumnGeneration(
            pricing_prob_solve_alg = Coluna.Algorithm.SolveIpForm(
                user_params = Coluna.Algorithm.UserOptimize(),
                moi_params = Coluna.Algorithm.MoiOptimize(
                    deactivate_artificial_vars = false,
                    enforce_integrality = false,
                ),
            ),
            stages_pricing_solver_ids = [1, 2, 3],
            throw_column_already_inserted_warning = true,
            smoothing_stabilization = 1.0,
        )
        colcutgen = Coluna.Algorithm.ColCutGenConquer(
            colgen = colgen,
            primal_heuristics = [],
            before_cutgen_user_algorithm = Coluna.Algorithm.BeforeCutGenAlgo(
                redcostfix_enum_algo, "Reduced cost fixing and enumeration",
            ),
            node_finalizer = Coluna.Algorithm.NodeFinalizer(
                solve_by_mip_algo, 0, "Solver by MIP",
            ),
            max_nb_cut_rounds = 9999,
        )
        branching = Coluna.Algorithm.StrongBranching()
        prodscore = Coluna.Algorithm.ProductScore()
        push!(
            branching.phases,
            Coluna.Algorithm.BranchingPhase(
                20, Coluna.Algorithm.RestrMasterLPConquer(), prodscore),
        )
        push!(branching.phases, Coluna.Algorithm.BranchingPhase(1, colcutgen, prodscore))
        push!(
            branching.rules,
            Coluna.Algorithm.Branching.PrioritisedBranchingRule(
                Coluna.Algorithm.SingleVarBranchingRule(), 1.0, 1.0),
        )

        tree_search = Coluna.Algorithm.TreeSearchAlgorithm(
            branchingtreefile = "BaPTree.dot",
            conqueralg = colcutgen,
            dividealg = branching,
        )
    end
    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(solver = tree_search),
        "default_optimizer" => () -> CPLEX.Optimizer(),
        # for the master & the subproblems
    )
    form = BlockModel(coluna) # , direct_model = true)
    if rcsp_path == ""
        push!(tree_search.opt, JuMP.unsafe_backend(form))
    end

    # Return the VrpSolver model containing the Coluna and RCSP models
    model = VrpModel(
        form, AffExpr(), RCSPProblem[],
        Vector{BlockDecomposition.Root{:VrpGraphs, Int64}}(undef, 1),
        Ptr{Cvoid}[], Ptr{Cvoid}(), 0, RCCPreSeparator[], CutCoeffManager(), VariableRef[],
        Dict{VariableRef, Int64}(), Dict{VariableRef, Vector{Int}}(), redcostfix_enum_algo, solve_by_mip_algo,
        VrpParameters[], Dict{String, Int}(), 0, "", Inf, Vector{Tuple{Int, Int}}[],
    )
    push!(tree_search.model_vec, model)
    return model
end

function set_branching_priority!(model::VrpModel, var_name::String, prior::Int)
    model.branch_priors[var_name] = prior
    return
end

function add_cut_callback!(::VrpModel, ::Function, ::String)
    @warn "add_cut_callback! is not implemented... ignoring it."
end
