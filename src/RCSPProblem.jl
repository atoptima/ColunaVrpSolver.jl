mutable struct RCSPProblem
    graph::VrpGraph
    solver::Ptr{Cvoid}
    buf_cint_1::Vector{Cint}
    buf_cint_2::Vector{Cint}
    buf_cint_3::Vector{Cint}
    buf_cint_4::Vector{Cint}
end

RCSPProblem(graph::VrpGraph) = RCSPProblem(
    graph, Ptr{Cvoid}(0), Cint[], Cint[], Cint[], Cint[]
)

struct PathVarData <: BlockDecomposition.AbstractCustomData
    graphid::Int
    arcids::Vector{Cint}
end

function build_solvers!(model::T) where {T <: AbstractVrpModel}
    # Instantiate an RCSP solver for each RCSP problem instance
    for prob in model.rcsp_instances
        prob.solver = ccall(
            (:createAndPrepareSolver_c, path), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}),
            prob.graph.cptr, get_rcsp_params(model.parameters[1], PARAM_CLASS_SOLVER)
        )
    end
    return
end

function run_rcsp_pricing(
    rcsp::RCSPProblem, phase::Int, var_rcosts::Vector{Float64}, r1cut_ptrs::Vector{Ptr{Cvoid}},
    r1cut_duals::Vector{Float64}
)
    # Call the RCSP pricing solver
    nb_r1cuts = length(r1cut_ptrs)
    output = [Ptr{Cvoid}(0)]
    nb_arcs = [Cint(0)]
    nb_paths = Int(ccall(
        (:runPricing_c, path), Cint,
        (
            Ptr{Cvoid}, Cint, Cint, Ptr{Float64}, Cint, Ref{Ptr{Cvoid}}, Ptr{Float64},
            Ref{Ptr{Cvoid}}, Ref{Cint}
        ),
        rcsp.solver, Cint(phase), Cint(length(var_rcosts)), var_rcosts, nb_r1cuts, r1cut_ptrs,
        r1cut_duals, output, nb_arcs
    ))

    # get the output paths, each path as a vertor of arc ids
    starts = rcsp.buf_cint_1
    arcs = rcsp.buf_cint_2
    resize!(starts, nb_paths + 1)
    resize!(arcs, Int(nb_arcs[1]))
    ccall(
        (:getOutputPaths_c, path), Cvoid, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}),
        output[1], starts, arcs
    )

    # convert and return them
    paths = Vector{Cint}[]
    for p in 1:nb_paths
        push!(paths, [a for a in arcs[(starts[p] + 1):starts[p + 1]]])
    end
    return paths
end

function run_rcsp_rcostfix_and_enum(
    rcsp::RCSPProblem, var_rcosts::Vector{Float64}, r1cut_ptrs::Vector{Ptr{Cvoid}},
    r1cut_duals::Vector{Float64}, threshold::Float64
)
    # Call the RCSP pricing solver to fix and enumerate
    nb_r1cuts = length(r1cut_ptrs)
    return (ccall(
        (:runRedCostFixingAndEnum_c, path), Cint,
        (Ptr{Cvoid}, Cint, Ptr{Float64}, Cint, Ref{Ptr{Cvoid}}, Ptr{Float64}, Float64),
        rcsp.solver, Cint(length(var_rcosts)), var_rcosts, nb_r1cuts, r1cut_ptrs, r1cut_duals,
        threshold
    ) != 0)
end

function record_rcsp_state(rcsp::RCSPProblem)
    # Call the RCSP pricing solver and return its current state
    return ccall(
        (:recordPricingState_c, path), Ptr{Cvoid}, (Ptr{Cvoid},), rcsp.solver
    )
end

function restore_rcsp_state(solver::Ptr{Cvoid}, state::Ptr{Cvoid})
    # Call the RCSP pricing solver to restore the state saved in `state`
    ccall(
        (:restorePricingState_c, path), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), solver, state
    )
    return
end

function release_rcsp_state(state::Ptr{Cvoid})
    # Call the RCSP pricing solver to release the memory used by the state saved at `state`
    ccall(
        (:releasePricingState_c, path), Cvoid, (Ptr{Cvoid},), state
    )
    return
end

function check_enumerated_paths(rcsp::RCSPProblem, paths::Vector{PathVarData})
    # Convet the paths data
    graphids = rcsp.buf_cint_1
    starts = rcsp.buf_cint_2
    arcids = rcsp.buf_cint_3
    empty!(graphids)
    empty!(starts)
    empty!(arcids)
    for p in paths
        push!(graphids, Cint(p.graphid))
        push!(starts, Cint(length(arcids)))
        for a in p.arcids
            push!(arcids, a)
        end
    end
    push!(starts, Cint(length(arcids)))

    # Call the RCSP pricing solver to select the relevant paths
    isrelevant = rcsp.buf_cint_4
    resize!(isrelevant, length(paths))
    ccall(
        (:checkEnumeratedPaths_c, path), Cvoid,
        (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        rcsp.solver, Cint(length(paths)), graphids, starts, arcids, isrelevant
    )
    return (isrelevant .!= 0)
end

function get_number_of_enum_paths(rcsp::RCSPProblem)
    return Int(ccall(
        (:getNumberOfEnumPaths_c, path), Cint, (Ptr{Cvoid},), rcsp.solver
    ))
end

function get_enumerated_paths(rcsp::RCSPProblem)
    # get the number of paths
    nb_paths = ccall(
        (:getNumberOfEnumPaths_c, path), Cint, (Ptr{Cvoid},), rcsp.solver
    )

    # get the total number of arcs in all paths
    starts = rcsp.buf_cint_1
    arcids = rcsp.buf_cint_2
    resize!(starts, nb_paths + 1)
    ccall(
        (:getEnumeratedPaths_c, path), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Cint, Ptr{Cint}),
        rcsp.solver, nb_paths, starts, 0, arcids
    )
    nb_arcs = starts[nb_paths + 1]

    # retrieve the paths
    resize!(arcids, nb_arcs)
    ccall(
        (:getEnumeratedPaths_c, path), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Cint, Ptr{Cint}),
        rcsp.solver, nb_paths, starts, nb_arcs, arcids
    )

    # convert and return
    return [[Int(a) for a in arcids[(starts[p] + 1):starts[p + 1]]] for p in 1:nb_paths]
end
