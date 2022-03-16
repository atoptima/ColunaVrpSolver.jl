module ColunaVrpSolver
    using Coluna, JuMP, BlockDecomposition, Gurobi

    export VrpModel, VrpGraph, add_resource!, add_graph!

    include("VrpModel.jl")
    include("RcspWrapper.jl")
end
