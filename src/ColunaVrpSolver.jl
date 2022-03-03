module ColunaVrpSolver
    using Coluna, JuMP, BlockDecomposition, Gurobi

    export VrpSolver

    include("Model.jl")
end
