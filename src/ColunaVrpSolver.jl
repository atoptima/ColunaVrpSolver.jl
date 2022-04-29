module ColunaVrpSolver
    using Coluna, JuMP, BlockDecomposition, Gurobi, MathOptInterface

    export VrpModel, VrpGraph, VrpOptimizer
    export add_resource!, set_resource_bounds!, add_arc!
    export set_arc_consumption!, add_graph!
    export set_cutoff!, optimize!

    include("VrpGraph.jl")
    include("VrpModel.jl")
    include("VrpOptimizer.jl")
end
