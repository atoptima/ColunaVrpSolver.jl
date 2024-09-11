module ColunaVrpSolver
using Coluna, JuMP, BlockDecomposition, CPLEX, MathOptInterface, Parameters, Printf

export VrpModel, VrpGraph, VrpOptimizer
export add_resource!, set_resource_bounds!
export add_arc!, add_arc_var_mapping!, set_arc_consumption!
export add_graph!, set_vertex_packing_sets!, define_elementarity_sets_distance_matrix!
export add_capacity_cut_separator!, set_branching_priority!, add_cut_callback!, set_cutoff!
export optimize!, get_objective_value, get_value

# Get the RCSP binary library path
global rcsp_path = get(ENV, "RCSP_LIB_PATH", "")

function dummy_ccall(_::Type{T}) where T
    if T == Cvoid
        return
    end
    return T(0)
end

macro try_ccall(func_and_path, rettype, args...)
    if rcsp_path == ""
        return esc(quote
            dummy_ccall($rettype)
        end)
    end
    args = map(esc, args)
    return quote
        ccall($func_and_path, $rettype, $(args...))
    end
end


include("Parameters.jl")
include("VrpGraph.jl")
include("RCSPProblem.jl")
include("VrpCuts.jl")
include("VrpModel.jl")
include("VrpMipSolver.jl")
include("VrpStorage.jl")
include("VrpOptimizer.jl")
end
