# Data structure to avoid cut coefficient repetitions minimizing allocations
mutable struct CutCoeffManager
    has_coeff::Vector{Int}
    nb_cuts::Int
end

CutCoeffManager() = CutCoeffManager(Int[], 0)

function nextcut!(ccm::CutCoeffManager, model::T) where {T<:AbstractVrpModel}
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
    bd_graphs::Vector{BlockDecomposition.Root{:VrpGraphs,Int64}}
    variables_by_id::Vector{VariableRef}
    branch_priors::Dict{String,Int}
    varids_by_var::Dict{VariableRef,Int}
    spids_by_var::Dict{VariableRef,Vector{Bool}}
    nb_subproblems::Int
    cfg_fname::String
    cutoffvalue::Float64
    packing_sets::Vector{Vector{Tuple{Int,Int}}}
    pset_to_id::Dict{Vector{Tuple{Int,Int}},Int}
    rcc_demands::Vector{Tuple{Cint,Vector{Cint}}}
    is_vertex_psets::Bool
end

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
    tree_search = BapcodTreeSearchWrapper(Coluna.Optimizer[], VrpModel[])
    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(solver=tree_search),
        "default_optimizer" => () -> CPLEX.Optimizer(),
        # for the master & the subproblems
    )
    form = BlockModel(coluna) # , direct_model = true)
    push!(tree_search.opt, JuMP.unsafe_backend(form))

    # Return the VrpSolver model containing the Coluna and RCSP models
    model = VrpModel(
        form, AffExpr(), RCSPProblem[],
        Vector{BlockDecomposition.Root{:VrpGraphs,Int64}}(undef, 1), VariableRef[], Dict{String,Int}(),
        Dict{VariableRef,Int64}(), Dict{VariableRef,Vector{Bool}}(), 0, "", Inf, Vector{Tuple{Int,Int}}[],
        Dict{Vector{Tuple{Int,Int}},Int}(), Tuple{Cint,Vector{Cint}}[], true,
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
