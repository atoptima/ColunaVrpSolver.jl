module ColunaVrpSolver
    using Coluna, JuMP, BlockDecomposition, Gurobi

    export VrpModel, VrpGraph, add_resource!, set_resource_bounds!, add_arc!
    export set_arc_consumption!, add_graph!

    include("RcspWrapper.jl")
    include("VrpModel.jl")
end
