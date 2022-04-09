module ColunaVrpSolver
    using Coluna, JuMP, BlockDecomposition, Gurobi, MathOptInterface

    export VrpModel, VrpGraph, add_resource!, set_resource_bounds!, add_arc!
    export set_arc_consumption!, add_graph!

    include("VrpGraph.jl")
    include("VrpModel.jl")
    include("VrpOptimizer.jl")
end
