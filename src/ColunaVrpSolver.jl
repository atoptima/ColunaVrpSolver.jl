__precompile__(false)
module ColunaVrpSolver
using Coluna, JuMP, BlockDecomposition, CPLEX, MathOptInterface, Parameters, Printf

export VrpModel, VrpGraph, VrpOptimizer
export add_resource!, set_resource_bounds!
export add_arc!, add_arc_var_mapping!, set_arc_consumption!
export add_graph!, set_vertex_packing_sets!, set_arc_packing_sets!, define_elementarity_sets_distance_matrix!
export add_elem_set_to_vertex_init_ng_neighbourhood!
export add_capacity_cut_separator!, set_branching_priority!, add_cut_callback!, set_cutoff!
export optimize!, get_objective_value, get_value

abstract type AbstractVrpModel end

include("Bapcod.jl")
include("VrpGraph.jl")
include("VrpModel.jl")
include("VrpOptimizer.jl")
end
