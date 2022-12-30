module ColunaVrpSolver
    using Coluna, JuMP, BlockDecomposition, Gurobi, MathOptInterface, Parameters, Printf

    export VrpModel, VrpGraph, VrpOptimizer
    export add_resource!, set_resource_bounds!
    export add_arc!, add_arc_var_mapping!, set_arc_consumption!
    export add_graph!, set_vertex_packing_sets!, define_elementarity_sets_distance_matrix!
    export add_capacity_cut_separator!, set_branching_priority!, add_cut_callback!, set_cutoff!
    export optimize!, get_objective_value, get_value

    include("Parameters.jl")
    include("VrpGraph.jl")
    include("RCSPProblem.jl")
    include("VrpCuts.jl")
    include("VrpModel.jl")
    include("VrpMipSolver.jl")
    include("VrpStorage.jl")
    include("VrpOptimizer.jl")
end
