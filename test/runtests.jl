using Test, ColunaVrpSolver, JuMP, Profile, PProf

include("integration/rcsp.jl")
include("integration/medsize.jl")

# run_rcsp_integration_tests()
run_mediumsize_instances_tests()
