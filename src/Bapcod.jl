bapcod_path = get(ENV, "BAPCOD_RCSP_LIB", "")

## ccall redef
macro bcm_ccall(func, args...)
    f = "bcInterfaceModel_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

macro bcr_ccall(func, args...)
    f = "bcRCSP_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

macro bcs_ccall(func, args...)
    f = "bcInterfaceSolve_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

macro bcsol_ccall(func, args...)
    f = "bcSolution_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

function new!(
    param_file::String,
    print_param::Bool,
    int_obj::Bool,
    int_valued_bound::Bool,
    argc::Cint,
    argv::Array{String, 1},
)
    @bcm_ccall("new", Ptr{Cvoid}, (Ptr{UInt8}, UInt8, UInt8, UInt8, Cint, Ptr{Ptr{UInt8}}),
        param_file, print_param, int_obj, int_valued_bound, argc, argv)
end

function init_model!(modelptr::Ptr{Cvoid}, nbrows::Cint, nbcols::Cint)
    @bcm_ccall("initModel", Cvoid, (Ptr{Cvoid}, Cint, Cint),
        modelptr, nbrows, nbcols)
end

function set_art_cost_value!(mptr::Ptr{Cvoid}, acv::Cdouble)
    @bcm_ccall("setArtCostValue", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, acv)
end

function set_obj_ub!(mptr::Ptr{Cvoid}, ub::Cdouble)
    @bcm_ccall("setObjUb", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, ub)
end

function register_sub_problem!(mptr::Ptr{Cvoid}, subproblemtype::Cint, spmid::Array)
    @bcm_ccall("registerSubProblem", Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}),
        mptr, subproblemtype, spmid)
end

function sptype_to_int(sp_type::Symbol)
    if sp_type == :MIP
        return 0
    elseif sp_type == :DW_MASTER || sp_type == :B_MASTER
        return 1
    elseif sp_type == :DW_SP
        return 2
    elseif sp_type == :B_SP
        return 3
    elseif sp_type == :ALL
        return 4
    else
        error(
            "Cannot recognize problem type : $sp_type. It must be :DW_MASTER or :DW_SP for Dantzig-Wolfe decomposition and :B_MASTER or :B_SP for Benders decomposition.",
        )
    end
end

function toArray(a)
    if isa(a, Integer)
        return [a]
    end
    arr = Vector{Int}()
    for i in a
        arr = vcat(arr, toArray(i))
    end
    return arr
end

function createMultiIndex(array_mid::Vector{Cint}, array::Vector{Int})
    length_array = length(array)
    length_array_mid = length(array_mid)
    (length(array) > 8) && error("BaPCod does not support multi-index with more than 8 indices.")
    for i in 1:length_array_mid
        if i <= length_array
            array_mid[i] = array[i]
        else
            array_mid[i] = (i == 1) ? 0 : -1
        end
    end
end

from_index_to_BaPCodindex(id, array_mid::Vector{Cint}) = createMultiIndex(array_mid, toArray(id))

function c_register_subproblems(mptr::Ptr{Cvoid}, spids)
    for (spid, sptype) in spids
        spmid = Array{Cint, 1}(undef, 8)
        from_index_to_BaPCodindex(spid, spmid)
        subproblemtype = Cint(sptype_to_int(sptype))
        register_sub_problem!(mptr, subproblemtype, spmid)
    end
end

function register_var!(
    mptr::Ptr{Cvoid},
    name::Symbol,
    column_id::Integer,
    sp_type::Integer,
    sp_mid::Array,
    var_mid::Array,
)
    @bcm_ccall("registerVar", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
        mptr, name, column_id, sp_type, sp_mid, var_mid)
end

function init_vars!(mptr::Ptr{Cvoid}, l::Array, u::Array, c::Array)
    @bcm_ccall("initVars", Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        mptr, l, u, c)
end

function c_register_vars(mptr::Ptr{Cvoid}, l, u, c, vars_decomposition)
    var_bcid = Array{Cint, 1}(undef, 8)
    sp_bcid = Array{Cint, 1}(undef, 8)

    for (column_id, (name, v_id, sp_type, sp_id)) in enumerate(vars_decomposition)
        # BaPCod needs an index
        from_index_to_BaPCodindex(v_id, var_bcid)
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        # Register the variable
        register_var!(mptr, name, Cint(column_id - 1), sptype_to_int(sp_type), sp_bcid, var_bcid)
    end
    init_vars!(mptr, l, u, c)
end

function init_cstrs!(mptr::Ptr{Cvoid}, starts::Array, rows_id::Array, nonzeros::Array, lb::Array, ub::Array)
    @bcm_ccall("initCstrs", Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        mptr, starts, rows_id, nonzeros, lb, ub)
end

function register_cstr!(
    mptr::Ptr{Cvoid},
    name::Symbol,
    row_id::Cint,
    sp_type::Integer,
    sp_mid::Array,
    cstr_mid::Array,
)
    @bcm_ccall("registerCstr", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
        mptr, string(name), row_id, sp_type, sp_mid, cstr_mid)
end

struct CMatrix
    starts::Array{Cint, 1}
    rows_id::Array{Cint, 1}
    nonzeros::Array{Cdouble, 1}
end

function c_register_cstrs(mptr::Ptr{Cvoid}, A, lb, ub, cstrs_decomposition)
    cstr_bcid = Array{Cint, 1}(undef, 8)
    sp_bcid = Array{Cint, 1}(undef, 8)

    for (row_id, (name, c_id, sp_type, sp_id)) in enumerate(cstrs_decomposition)
        # BaPCod needs an idnex
        from_index_to_BaPCodindex(c_id, cstr_bcid)
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        # Register the constraint
        register_cstr!(mptr, name, Cint(row_id - 1), sptype_to_int(sp_type), sp_bcid, cstr_bcid)
    end
    init_cstrs!(mptr, A.starts, A.rows_id, A.nonzeros, lb, ub)
end

function sub_problem_mult!(mptr::Ptr{Cvoid}, mult_lb::Cint, mult_ub::Cint, sp_type::Integer, sp_bcid::Vector{Cint})
    @bcm_ccall("subProblemMult", Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cint}),
        mptr, mult_lb, mult_ub, sp_type, sp_bcid)
end

function c_set_sp_multiplicities(mptr::Ptr{Cvoid}, sp_mult)
    sp_bcid = Array{Cint, 1}(undef, 8)
    for (sp_id, sp_type, mult_lb, mult_ub) in sp_mult
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        status = sub_problem_mult!(mptr, Cint(mult_lb), Cint(mult_ub), sptype_to_int(sp_type), sp_bcid)
        (status != 1) &&
            error("Cannot set multiplicity on the subproblem with the index $sp_id. Make sure it exists.")
    end
end

function set_var_priority_in_master!(
    modelptr::Ptr{Cvoid},
    varname::Symbol,
    sp_bctype::Integer,
    sp_bcid::Array,
    priority::Cdouble,
)
    @bcm_ccall("setVarPriorityInMaster", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Cdouble),
        modelptr, varname, sp_bctype, sp_bcid, priority)
end

function c_vars_branching_priorities(modelptr::Ptr{Cvoid}, p::Vector{Tuple{Symbol, Symbol, Int, Cdouble}})
    sp_bcid = Array{Cint, 1}(undef, 8)
    for (varname, sp_type, sp_id, priority) in p
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        sp_bctype = sptype_to_int(sp_type)
        status = set_var_priority_in_master!(modelptr, varname, sp_bctype, sp_bcid, priority)
        (status != 1) && error("Cannot set branching priority on variables named $varname.")
    end
end

function wbcr_new(
    c_model::Ptr{Cvoid},
    sp_bctype::Integer,
    sp_bcid::Array,
    nb_nodes::Integer,
    nb_es::Integer,
    nb_ps::Integer,
    nb_cs::Integer,
)
    ptr = @bcr_ccall("new", Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Ptr{Int}, Cint, Cint, Cint, Cint),
        c_model, Cint(sp_bctype), sp_bcid, Cint(nb_nodes), Cint(nb_es), Cint(nb_ps), Cint(nb_cs))
end

function new_network!(
    c_model::Ptr{Cvoid},
    sp_id::Int,
    sp_type::Symbol,
    nb_nodes::Int,
    nb_es::Int,
    nb_ps::Int,
    nb_cs::Int,
)
    sp_bcid = Array{Cint}(undef, 8)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    sp_bctype = sptype_to_int(sp_type)
    return wbcr_new(c_model, sp_bctype, sp_bcid, nb_nodes, nb_es, nb_ps, nb_cs)
end

function wbcr_new_resource(c_net::Ptr{Cvoid}, res_id::Integer)
    status = @bcr_ccall("newResource", Cint, (Ptr{Cvoid}, Cint),
        c_net, Cint(res_id))
    (status != 1) && error("Cannot create resource $res_id.")
end

function wbcr_set_vertex_consumption_lb(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, lb::Cdouble)
    @bcr_ccall("setVertexConsumptionLB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(n_id), Cint(res_id), lb)
end

function wbcr_set_vertex_consumption_ub(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, ub::Cdouble)
    @bcr_ccall("setVertexConsumptionUB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(n_id), Cint(res_id), ub)
end

function wbcr_attach_elementarity_set_to_node(c_net::Ptr{Cvoid}, n_id::Integer, es_id::Integer)
    @bcr_ccall("attachElementaritySetToNode", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(n_id), Cint(es_id))
end

function wbcr_add_vertex_to_mem_of_elementarity_set(c_net::Ptr{Cvoid}, n_id::Integer, es_id::Integer)
    status = @bcr_ccall("addVertexToMemOfElementaritySet", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(n_id), Cint(es_id))
    (status != 1) && error("Cannot add vertex $n_id to memory of elementarity set $n_es.")
end

function wbcr_add_vertex_to_packing_set(c_net::Ptr{Cvoid}, n_id::Integer, ps_id::Integer)
    status = @bcr_ccall("addVertexToPackingSet", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(n_id), Cint(ps_id))
end

function wbc_add_generic_lim_mem_one_cut(c_model::Ptr{Cvoid})
    status = @bcr_ccall("addGenericLimMemOneCut", Cint, (Ptr{Cvoid},),
        c_model)
    (status != 1) && error("Cannot add the generic lim-mem-one-rank cut.")
end

function wbcr_set_source(c_net::Ptr{Cvoid}, n_id::Integer)
    @bcr_ccall("setSource", Cint, (Ptr{Cvoid}, Cint), c_net, Cint(n_id))
end

function wbcr_set_sink(c_net::Ptr{Cvoid}, n_id::Integer)
    @bcr_ccall("setSink", Cint, (Ptr{Cvoid}, Cint), c_net, Cint(n_id))
end

function wbcr_new_arc(c_net::Ptr{Cvoid}, src::Integer, dst::Integer, cost::Cdouble)
    edge_id = @bcr_ccall("newArc", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(src), Cint(dst), cost)
    return edge_id
end

function wbcr_attach_bcvar_to_arc(c_net::Ptr{Cvoid}, edge_id::Integer, c_model::Ptr{Cvoid}, var_col::Cint)
    @bcr_ccall("attachBcVarToArc", Cint, (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Cint, Cdouble),
        c_net, Cint(edge_id), c_model, var_col, Cdouble(1.0))
end

function wbcr_set_edge_consumption_value(c_net::Ptr{Cvoid}, edge_id::Integer, res_id::Integer, value::Cdouble)
    @bcr_ccall("setEdgeConsumptionValue", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(edge_id), Cint(res_id), value)
end

function wbcr_create_oracle(
    c_net::Ptr{Cvoid},
    c_model::Ptr{Cvoid},
    sp_type::Integer,
    sp_id::Array,
    save_standalone::Bool,
    standalone_filename::String,
)
    @bcr_ccall("createOracle", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cint}, UInt8, Ptr{UInt8}),
        c_net, c_model, sp_type, sp_id, save_standalone, standalone_filename)
end

function new_oracle!(c_net::Ptr{Cvoid}, c_model::Ptr{Cvoid}, sp_type::Symbol, sp_id::Int)
    sp_bcid = Array{Cint, 1}(undef, 8)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    sp_bctype = sptype_to_int(sp_type)
    wbcr_create_oracle(c_net, c_model, sp_bctype, sp_bcid, false, "")
end

function c_optimize(modelptr::Ptr{Cvoid}, solution::Ptr{Cvoid})
    @bcs_ccall("optimize", Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), modelptr, solution)
end

function new_sol!()
    @bcsol_ccall("new", Ptr{Cvoid}, ())
end

function c_getValues(mptr::Ptr{Cvoid}, solution::Ptr{Cvoid}, vector, nbvars::Integer)
    @bcsol_ccall("getValues", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cdouble}, Cint),
        mptr, solution, vector, nbvars)
end

function c_getMultiplicity(solution::Ptr{Cvoid}, mult)
    @bcsol_ccall("getMultiplicity", Cint, (Ptr{Cvoid}, Ref{Cint}),
        solution, mult)
end

function c_start(solution::Ptr{Cvoid}, mptr)
    @bcsol_ccall("start", Cint, (Ptr{Cvoid}, Ptr{Cvoid}), solution, mptr)
end

function c_next(solution::Ptr{Cvoid})
    @bcsol_ccall("next", Cint, (Ptr{Cvoid},), solution)
end

# Build the solution stored in a matrix. Each row is a column.
function get_solution(modelptr::Ptr{Cvoid}, solptr::Ptr{Cvoid}, nbvars::Int)
    # For each variable we create an array containing the variable value for each solution
    sol = Tuple{Int, Array{Tuple{Int, Float64}}}[]
    status = c_start(solptr, modelptr)
    while status == 1
        vector = Array{Cdouble}(undef, nbvars)
        c_getValues(modelptr, solptr, vector, nbvars)
        mult = Ref{Cint}(0)
        nonzeros = [(i, vector[i]) for i in 1:nbvars if vector[i] != 0.0]
        if !isempty(nonzeros)
            status = c_getMultiplicity(solptr, mult)
            push!(sol, (Int(mult[]), nonzeros))
        end
        status = c_next(solptr)
    end
    return sol
end

function convert_solution(
    masterform::Coluna.MathProg.Formulation,
    sol::Vector{Tuple{Int, Array{Tuple{Int, Float64}}}},
    colid_to_varid::Vector{Coluna.MathProg.Id{Coluna.MathProg.Variable}},
    colid_to_cost::Vector{Float64},
    colid_to_spform::Vector{Vector{Coluna.MathProg.Formulation{Coluna.MathProg.DwSp}}},
)
    varids = Coluna.MathProg.VarId[]
    varcoeffs = Float64[]
    cost = 0.0
    for (mult, subsol) in sol
        spform = colid_to_spform[subsol[1][1]][1]

        # build the subproblem solution data
        subvarids = [colid_to_varid[colid] for (colid, _) in subsol]
        subvarcoeffs = [coeff for (_, coeff) in subsol]
        push!(subvarids, spform.duty_data.setup_var)
        push!(subvarcoeffs, 1.0)
        subcost = sum(coeff * colid_to_cost[colid] for (colid, coeff) in subsol)

        # add the subproblem solution to the subproblem
        subsol = Coluna.MathProg.PrimalSolution(
            spform, subvarids, subvarcoeffs, subcost, Coluna.FEASIBLE_SOL,
        )
        col_id = Coluna.MathProg.insert_column!(masterform, subsol, "MC")
        mc_var = Coluna.MathProg.getvar(masterform, col_id)

        # accumulate the master solution data
        push!(varids, Coluna.MathProg.getid(mc_var))
        push!(varcoeffs, Float64(mult))
        cost += Float64(mult) * subcost
    end

    # add the master solution to the master problem and return it
    return Coluna.MathProg.PrimalSolution(
        masterform, varids, varcoeffs, cost, Coluna.FEASIBLE_SOL,
    )
end

struct BapcodTreeSearchWrapper{M <: AbstractVrpModel} <: Coluna.Algorithm.AbstractOptimizationAlgorithm
    opt::Vector{Coluna.Optimizer}
    model_vec::Vector{M}
end

function Coluna.Algorithm.run!(
    algo::BapcodTreeSearchWrapper,
    env::Coluna.Env,
    reform::Coluna.MathProg.Reformulation,
    input::Coluna.Algorithm.OptimizationState,
)
    # Get the subproblems sorted by id (FIXME: will not work if the user creates other subproblems)
    sps = [form for (_, form) in Coluna.MathProg.get_dw_pricing_sps(reform)]
    sort!(sps, by = form -> Coluna.ColunaBase.getuid(form))

    model = algo.model_vec[1]
    masterform = reform.master
    # print("$(Coluna.MathProg.getobjsense(masterform)) ")
    first = true
    varid_to_prior = Dict{Coluna.MathProg.Id{Coluna.MathProg.Variable}, Float64}()
    varid_to_varref = Dict{Coluna.MathProg.Id{Coluna.MathProg.Variable}, VariableRef}()
    for varref in JuMP.all_variables(model.formulation)
        vid = Coluna._get_varid_of_origvar_in_form(algo.opt[1].env, masterform, JuMP.index(varref))
        varid_to_prior[vid] = 0.0
        varid_to_varref[vid] = varref
    end
    for (varname, prior) in model.branch_priors
        for varref in model.formulation[Symbol(varname)]
            vid = Coluna._get_varid_of_origvar_in_form(algo.opt[1].env, masterform, JuMP.index(varref))
            varid_to_prior[vid] = Float64(prior)
        end
    end
    ncols = Cint(0) # number of columns in the block-diagonal matrix of the original problem
    lbs = Cdouble[]
    ubs = Cdouble[]
    costs = Cdouble[]
    vars = Tuple{Symbol, Int, Symbol, Int}[]
    priors = Tuple{Symbol, Symbol, Int, Cdouble}[]
    varid_to_colids = Dict{Coluna.MathProg.Id{Coluna.MathProg.Variable}, Vector{Cint}}()
    colid_to_varid = Coluna.MathProg.Id{Coluna.MathProg.Variable}[]
    colid_to_spform = Vector{Coluna.MathProg.Formulation{Coluna.MathProg.DwSp}}[]
    for (var_id, var) in Coluna.MathProg.getvars(masterform)
        if Coluna.MathProg.getduty(var_id) <= Coluna.MathProg.MasterPureVar ||
           Coluna.MathProg.getduty(var_id) <= Coluna.MathProg.MasterRepPricingVar
            name = Coluna.MathProg.getname(masterform, var_id)
            if var.curdata.cost != 0
                if first
                    # print("$(var.curdata.cost) * $(name)")
                    first = false
                else
                    # print(" + $(var.curdata.cost) * $(name)")
                end
            end
            spids = get(model.spids_by_var, varid_to_varref[var_id], Bool[])
            if isempty(spids)
                push!(lbs, Cdouble(var.curdata.lb))
                push!(ubs, Cdouble(var.curdata.ub))
                push!(costs, Cdouble(var.curdata.cost))
                push!(vars, (Symbol(name), Int(ncols), :DW_MASTER, 0))
                push!(priors, (Symbol(name), :DW_MASTER, 0, varid_to_prior[var_id]))
                push!(colid_to_varid, var_id)
                push!(colid_to_spform, Coluna.MathProg.Id{Coluna.MathProg.Variable}[])
                ncols += Cint(1)
            else
                for (spid1, used) in enumerate(spids)
                    !used && continue
                    spid = spid1 - 1
                    push!(lbs, Cdouble(0.0)) # var.curdata.lb))
                    push!(ubs, Cdouble(Inf)) # var.curdata.ub))
                    push!(costs, Cdouble(var.curdata.cost))
                    varsymbol = Symbol(name * "_$spid")
                    push!(vars, (varsymbol, Int(ncols), :DW_SP, spid))
                    push!(priors, (varsymbol, :DW_SP, spid, varid_to_prior[var_id]))
                    if haskey(varid_to_colids, var_id)
                        push!(varid_to_colids[var_id], ncols)
                    else
                        varid_to_colids[var_id] = [ncols]
                    end
                    push!(colid_to_varid, var_id)
                    push!(colid_to_spform, [sps[spid1]])
                    ncols += Cint(1)
                end
            end
        end
    end
    # println("\n")
    # @show vars
    # @show lbs
    # @show ubs
    matrix = Coluna.MathProg.getcoefmatrix(masterform)
    nconstrs = Cint(0)
    clbs = Cdouble[]
    cubs = Cdouble[]
    constrs = Tuple{Symbol, Int, Symbol, Int}[]
    constr_id_to_row_id = Dict{Coluna.MathProg.Id{Coluna.MathProg.Constraint}, Cint}()
    for (constr_id, constr) in Coluna.MathProg.getconstrs(masterform)
        if Coluna.MathProg.getduty(constr_id) <= Coluna.MathProg.AbstractMasterOriginConstr
            name = Coluna.MathProg.getname(masterform, constr_id)
            # print("$(name) ($(Coluna.MathProg.getduty(constr_id))):")
            # first = true
            # for (var_id, coeff) in @view matrix[constr_id, :]
            #     varname = Coluna.MathProg.getname(masterform, var_id)
            #     if Coluna.MathProg.getduty(var_id) <= Coluna.MathProg.MasterPureVar ||
            #        Coluna.MathProg.getduty(var_id) <= Coluna.MathProg.MasterRepPricingVar
            #         if first
            #             print(" $coeff * $varname")
            #             first = false
            #         else
            #             print(" + $coeff * $varname")
            #         end
            #     end
            # end
            # sense = constr.curdata.sense
            # rhs = constr.curdata.rhs
            # if sense == Coluna.MathProg.Less
            #     println(" <= $rhs")
            # elseif sense == Coluna.MathProg.Greater
            #     println(" >= $rhs")
            # else
            #     println(" = $rhs")
            # end
            if constr.curdata.sense == Coluna.MathProg.Less
                push!(clbs, Cdouble(-Inf))
                push!(cubs, Cdouble(constr.curdata.rhs))
            elseif constr.curdata.sense == Coluna.MathProg.Greater
                push!(clbs, Cdouble(constr.curdata.rhs))
                push!(cubs, Cdouble(Inf))
            elseif constr.curdata.sense == Coluna.MathProg.Equal
                push!(clbs, Cdouble(constr.curdata.rhs))
                push!(cubs, Cdouble(constr.curdata.rhs))
            else
                error("Cannot recognize constraint sense.")
            end
            push!(constrs, (Symbol(name), Int(nconstrs), :DW_MASTER, -1))
            constr_id_to_row_id[constr_id] = nconstrs
            nconstrs += Cint(1)
        end
    end
    starts = Cint[]
    rows_id = Cint[]
    nonzeros = Cdouble[]
    for (var_id, _) in Coluna.MathProg.getvars(masterform)
        if Coluna.MathProg.getduty(var_id) <= Coluna.MathProg.MasterPureVar ||
           Coluna.MathProg.getduty(var_id) <= Coluna.MathProg.MasterRepPricingVar
            nb_var_cols = length(get(model.spids_by_var, varid_to_varref[var_id], [0]))
            for _ in 1:nb_var_cols
                push!(starts, Cint(length(nonzeros)))
                for (constr_id, coeff) in @view matrix[:, var_id]
                    if Coluna.MathProg.getduty(constr_id) <= Coluna.MathProg.AbstractMasterOriginConstr
                        push!(rows_id, constr_id_to_row_id[constr_id])
                        push!(nonzeros, Cdouble(coeff))
                    end
                end
            end
        end
    end
    push!(starts, Cint(length(nonzeros)))
    # @show starts
    # @show rows_id
    # @show nonzeros
    model_ptr = new!(model.cfg_fname, true, true, false, Cint(0), String[])
    init_model!(model_ptr, nconstrs, ncols)
    set_art_cost_value!(model_ptr, Cdouble(10000))
    set_obj_ub!(model_ptr, Cdouble(model.cutoffvalue))
    c_register_subproblems(model_ptr, [(spid, :DW_SP) for spid in 0:(length(model.rcsp_instances)-1)])
    c_register_vars(model_ptr, lbs, ubs, costs, vars)
    c_register_cstrs(model_ptr, CMatrix(starts, rows_id, nonzeros), clbs, cubs, constrs)
    c_set_sp_multiplicities(
        model_ptr,
        [
            (spid, :DW_SP, model.rcsp_instances[spid+1].graph.bounds...) for
            spid in 0:(length(model.rcsp_instances)-1)
        ],
    )
    c_vars_branching_priorities(model_ptr, priors)

    # Build the RCSP network for each subproblem
    net_ptrs = Ptr{Cvoid}[]
    for spid in 0:(length(model.rcsp_instances)-1)
        graph = model.rcsp_instances[spid+1].graph
        nb_nodes = maximum(graph.vert_ids) + 1
        nb_psets = length(model.packing_sets)
        nb_elemsets = length(graph.elem_sets)
        c_net_ptr = new_network!(model_ptr, spid, :DW_SP, nb_nodes, nb_psets, nb_elemsets, 0)
        push!(net_ptrs, c_net_ptr)
        for resid in 0:(graph.nb_resources-1)
            wbcr_new_resource(c_net_ptr, resid)
            for i in 1:nb_nodes
                i_ = (i == nb_nodes) ? 1 : i
                wbcr_set_vertex_consumption_lb(
                    c_net_ptr,
                    i - 1,
                    resid,
                    Cdouble(graph.res_bounds[i_][resid+1][1]),
                )
                wbcr_set_vertex_consumption_ub(
                    c_net_ptr,
                    i - 1,
                    resid,
                    Cdouble(graph.res_bounds[i_][resid+1][2]),
                )
                # @show (i - 1), graph.res_bounds[i_][resid+1]
            end
        end
        for es_id in eachindex(graph.elem_sets)
            elem_set = graph.elem_sets[es_id]
            for i in elem_set
                j = graph.vert_ids[i+1]
                wbcr_attach_elementarity_set_to_node(c_net_ptr, j, es_id - 1)
            end
            dists = graph.dist_matrix[es_id]
            neighs = [k for k in 0:(nb_nodes-1) if k != graph.src_id && k != graph.snk_id]
            sort!(neighs, by = x -> dists[x])
            for (k, j) in enumerate(neighs)
                wbcr_add_vertex_to_mem_of_elementarity_set(c_net_ptr, j, es_id - 1)
                if k == model.parameters[1].coluna_vrp_params.RCSPmaxNGneighbourhoodSize
                    break
                end
            end
        end
        wbcr_set_source(c_net_ptr, graph.src_id)
        wbcr_set_sink(c_net_ptr, graph.snk_id)
        # println("Mappings and Consumptions:")
        for (id1, (tail, head)) in enumerate(graph.arcs)
            push!(
                graph.arc_ids,
                wbcr_new_arc(c_net_ptr, graph.vert_ids[tail+1], graph.vert_ids[head+1], Cdouble(0.0)),
            )
            for var in graph.mappings[id1]
                vid = Coluna._get_varid_of_origvar_in_form(algo.opt[1].env, masterform, JuMP.index(var))
                colids = varid_to_colids[vid]
                for colid in colids
                    # print(" $(Coluna.MathProg.getname(masterform, vid)), $colid -> $((tail, head)),")
                    wbcr_attach_bcvar_to_arc(c_net_ptr, graph.arc_ids[id1], model_ptr, colid)
                end
            end
            for resid in 0:(graph.nb_resources-1)
                wbcr_set_edge_consumption_value(
                    c_net_ptr,
                    graph.arc_ids[id1],
                    resid,
                    Cdouble(graph.res_cons[id1][resid+1]),
                )
                # print(" $(graph.res_cons[id1][resid+1]),")
            end
        end
        # println("\n")
        new_oracle!(c_net_ptr, model_ptr, :DW_SP, spid)
    end

    # Set the packing sets
    for (ps_id1, ps) in enumerate(model.packing_sets)
        for (spid, i) in ps
            graph = model.rcsp_instances[spid+1].graph
            wbcr_add_vertex_to_packing_set(net_ptrs[spid+1], graph.vert_ids[i+1], ps_id1 - 1)
        end
    end
    if !isempty(model.packing_sets)
        wbc_add_generic_lim_mem_one_cut(model_ptr)
    end

    sol_ptr = new_sol!()
    c_optimize(model_ptr, sol_ptr)
    sol = get_solution(model_ptr, sol_ptr, Int(ncols))
    # println("Solution: ", sol)

    # Set the optimization state to output
    output = Coluna.Algorithm.OptimizationState(
        masterform,
        ip_primal_bound = Coluna.Algorithm.get_ip_primal_bound(input),
        termination_status = Coluna.OPTIMAL,
    )

    # Build the primal solution if any
    if !isempty(sol)
        Coluna.Algorithm.add_ip_primal_sol!(
            output,
            convert_solution(masterform, sol, colid_to_varid, costs, colid_to_spform),
        )
        optcost = Coluna.ColunaBase.getvalue(Coluna.Algorithm.get_ip_primal_bound(output))
        Coluna.Algorithm.set_ip_dual_bound!(output, Coluna.MathProg.DualBound(masterform, optcost))
    end

    # Return the optimization state
    return output
end
