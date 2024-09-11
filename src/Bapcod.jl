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
