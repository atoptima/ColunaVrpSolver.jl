using ColunaVrpSolver

@testset "Calling the C++ library to solve the RCSP problem" begin
    @test haskey(ENV, "RCSP_LIB_PATH")
    VrpSolver()
end
