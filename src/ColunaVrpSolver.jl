module ColunaVrpSolver
    using Coluna, JuMP, BlockDecomposition, Gurobi, MathOptInterface

    export VrpModel, VrpGraph, VrpOptimizer
    export add_resource!, set_resource_bounds!
    export add_arc!, add_arc_var_mapping!, set_arc_consumption!
    export add_graph!, set_cutoff!, optimize!
    export get_objective_value, get_value

    include("VrpGraph.jl")
    include("VrpModel.jl")
    include("VrpOptimizer.jl")
end
