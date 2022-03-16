# Get the binary library path
global path = "$(ENV["RCSP_LIB_PATH"])"

mutable struct VrpGraph
    cptr::Ptr{Cvoid}
    bounds::Tuple{Float64, Float64}
    next_resid::Int
end

function VrpGraph(model::VrpModel, vertices::Vector{Int}, source::Int, sink::Int, bounds::Tuple{Int, Int})
    cptr_ = ccall(
        (:createGraph_c, path), Ptr{Cvoid}, (Cint, Ptr{Cint}, Cint, Cint),
        Cint(length(vertices)), [Cint(v) for v in vertices], Cint(source), Cint(sink)
    )
    graph = VrpGraph(cptr_, Float64.(bounds), 0)
    return graph
end

function add_resource!(graph::VrpGraph; main = false)
    id = graph.next_resid
    graph.next_resid += 1
    ccall((:addResource_c, path), Cvoid, (Ptr{Cvoid}, Cint, Cint), graph.cptr, Cint(id), Cint(main))
end

function add_graph!(model::VrpModel, graph::VrpGraph)
    # Instantiate an RCSP solver
    rcsp = ccall((:createAndPrepareSolver_c, path), Ptr{Cvoid}, (Ptr{Cvoid},), graph.cptr)

    # Add it to the VRP model
    push!(model.rcsp_instances, rcsp)
end

