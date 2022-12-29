@with_kw mutable struct ColunaVrpParams
    StrongBranchingPhaseOneCandidatesNumber::Int = 100
    StrongBranchingPhaseOneTreeSizeEstimRatio::Float64 = 0.3
    StrongBranchingPhaseTwoCandidatesNumber::Int = 3
    StrongBranchingPhaseTwoTreeSizeEstimRatio::Float64 = 0.1
    CutTailingOffThreshold::Float64 = 0.02
    ReducedCostFixingThreshold::Float64 = 0.9
    RCSPmaxNumOfEnumSolutionsForMIP::Int = 10000
    relOptimalityGapTolerance::Float64 = 1e-6
end

function print_params(params::ColunaVrpParams)
    println("")
    println("Values of ColunaVrpSolverParameters:")
    println("====================================================")
    for param in fieldnames(ColunaVrpParams)
        println(param, " = ", getfield(params, param))
    end
    return nothing
end

mutable struct VrpParameters
    rcsp_params::Vector{Ptr{Cvoid}}
    coluna_vrp_params::ColunaVrpParams
end

const PARAM_CLASS_COLUNA::Cint = -1
const PARAM_CLASS_ROUND_CAP_CUTS_SEPARATOR::Cint = 0
const PARAM_CLASS_STRONG_KPATH_SEPARATOR::Cint = 1
const PARAM_CLASS_ROUTE_LOAD_KNAPSACK_CUTS_SEPARATOR::Cint = 2
const PARAM_CLASS_LIM_MEM_RANK_ONE_CUTS_SEPARATOR::Cint = 3
const PARAM_CLASS_SOLVER::Cint = 4
const MAX_PARAM_CLASSES::Cint = 5

get_rcsp_params(params::Vector{Ptr{Cvoid}}, classid::Cint) = params[classid + 1]
get_rcsp_params(params::VrpParameters, classid::Cint) =
    get_rcsp_params(params.rcsp_params, classid)

# ====== This code has been copied from `https://rosettacode.org/`
function striplinecomment(a::String, cchars::String="#;")
    b = strip(a)
    0 < length(cchars) || return b
    for c in cchars
        r = Regex(@sprintf "\\%c.*" c)
        b = replace(b, r => "")
    end
    strip(b)
end
# ======

is_coluna_param(name::AbstractString) = Symbol(name) in fieldnames(ColunaVrpParams)
is_bool(val::AbstractString) = (String(val) in ["false", "true"])
is_int(val::AbstractString) = all(isnothing.(findfirst.(['.', 'e', 'E'], val)))

function setparam!(
    params_class::Cint, coluna_vrp_params::ColunaVrpParams, rcsp_params::Vector{Ptr{Cvoid}},
    param_name::String, value::T
) where T
    if params_class == PARAM_CLASS_COLUNA
        p = Symbol(param_name)
        if p in fieldnames(ColunaVrpParams)
            if typeof(getfield(coluna_vrp_params, p)) == T
                setfield!(coluna_vrp_params, p, value)
            else
                @warn "Invalid $T parameter value for $param_name in file $fname"
            end
        else
            @warn "Unknown parameter $param_name in file $fname"
        end
    else
        if !set_rcsp_parameter(get_rcsp_params(rcsp_params, params_class), param_name, value)
            @warn "Unknown parameter $param_name in file $fname"
        end
    end
    return
end

function VrpParameters(fname::String)
    rcsp_params = [
        create_rcsp_parameters(params_class)
        for params_class in Cint(0):(MAX_PARAM_CLASSES - Cint(1))
    ]
    coluna_vrp_params = ColunaVrpParams()
    buf = split(striplinecomment(read(fname, String)), ('=', ' ', '\n'), keepempty = false)
    for line in 1:2:length(buf)
        if buf[line][1:4] == "RCSP" && !is_coluna_param(buf[line])
            if buf[line][5:11] == "rankOneCuts"
                params_class = PARAM_CLASS_LIM_MEM_RANK_ONE_CUTS_SEPARATOR
                param_name = lowercasefirst(String(buf[line][12:end]))
            else
                params_class = PARAM_CLASS_SOLVER
                param_name = String(buf[line][5:end])
            end
        else
            params_class = PARAM_CLASS_COLUNA
            param_name = String(buf[line])
        end
        if is_bool(buf[line + 1])
            setparam!(
                params_class, coluna_vrp_params, rcsp_params, param_name,
                parse(Bool, buf[line + 1])
            )
        elseif is_int(buf[line + 1])
            setparam!(
                params_class, coluna_vrp_params, rcsp_params, param_name,
                parse(Int, buf[line + 1])
            )
        else
            setparam!(
                params_class, coluna_vrp_params, rcsp_params, param_name,
                parse(Float64, buf[line + 1])
            )
        end
    end
    return VrpParameters(rcsp_params, coluna_vrp_params)
end

function create_rcsp_parameters(class::Cint)
    # Create an instance of a parameters class
    return ccall(
        (:createParameters_c, path), Ptr{Cvoid}, (Cint,), class
    )
end

function set_rcsp_parameter(params::Ptr{Cvoid}, name::String, value::Bool)
    # Set the parameter with name `name` to the value `value` in `params`
    return (ccall(
        (:setBoolParamValue_c, path), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring, Cint),
        params, name, value ? Cint(1) : Cint(0)
    ) != 0)
end

function set_rcsp_parameter(params::Ptr{Cvoid}, name::String, value::Int)
    # Set the parameter with name `name` to the value `value` in `params`
    return (ccall(
        (:setIntParamValue_c, path), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring, Cint),
        params, name, Cint(value)
    ) != 0)
end

function set_rcsp_parameter(params::Ptr{Cvoid}, name::String, value::Float64)
    # Set the parameter with name `name` to the value `value` in `params`
    return (ccall(
        (:setDoubleParamValue_c, path), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring, Float64),
        params, name, value
    ) != 0)
end
