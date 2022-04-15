# Get the binary library path
global path = "$(ENV["RCSP_LIB_PATH"])"

abstract type AbstractVrpModel end

mutable struct VrpGraph
    id::Int
    cptr::Ptr{Cvoid}
    bounds::Tuple{Float64, Float64}
    orig_sink::Int
    new_sink::Int
    mappings::Dict{Int, Vector{VariableRef}}
    max_arcid::Int
end

function VrpGraph(
    _::T, vertices::Vector{Int}, source::Int, sink::Int, bounds::Tuple{Int, Int}
) where T
    new_sink = sink
    vertices_ = copy(vertices)
    if sink == source
        new_sink = maximum(vertices) + 1
        push!(vertices_, new_sink)
    end
    cptr_ = ccall(
        (:createGraph_c, path), Ptr{Cvoid}, (Cint, Ptr{Cint}, Cint, Cint),
        Cint(length(vertices_)), [Cint(v) for v in vertices_], Cint(source), Cint(new_sink)
    )
    graph = VrpGraph(
        0, cptr_, Float64.(bounds), sink, new_sink, Dict{Int, Vector{VariableRef}}(), 0
    )
    return graph
end

function add_resource!(graph::VrpGraph; main = false)
    id = ccall((:addResource_c, path), Cint, (Ptr{Cvoid}, Cint), graph.cptr, Cint(main))
    return Int(id)
end

function set_resource_bounds!(graph::VrpGraph, vertid::Int, resid::Int, lb::Float64, ub::Float64)
    vid = (vertid == graph.orig_sink) ? graph.new_sink : vertid
    ccall(
        (:setResourceBounds_c, path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64, Float64),
        graph.cptr, Cint(vid), Cint(resid), lb, ub
    )
    return
end

function add_arc!(graph::VrpGraph, tail::Int, head::Int)
    h = (head == graph.orig_sink) ? graph.new_sink : head
    id = ccall(
        (:addArc_c, path), Cint, (Ptr{Cvoid}, Cint, Cint), graph.cptr, Cint(tail), Cint(h)
    )
    graph.max_arcid = max(graph.max_arcid, Int(id))
    return Int(id)
end

function set_arc_consumption!(graph::VrpGraph, arcid::Int, resid::Int, cons::Float64)
    ccall(
        (:setArcConsumption_c, path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64),
        graph.cptr, Cint(arcid), Cint(resid), cons
    )
    return
end

function add_arc_var_mapping!(graph::VrpGraph, arcid::Int, var::VariableRef)
    mapped = get(graph.mappings, arcid, VariableRef[])
    if isempty(mapped)
        graph.mappings[arcid] = mapped
    end
    push!(mapped, var)
end

function add_graph!(model::T, graph::VrpGraph) where {T <: AbstractVrpModel}
    # Add the graph to the VRP model
    push!(model.graphs, graph)
    graph.id = length(model.graphs)

    # Instantiate an RCSP solver
    rcsp = ccall((:createAndPrepareSolver_c, path), Ptr{Cvoid}, (Ptr{Cvoid},), graph.cptr)

    # Add it to the VRP model
    push!(model.rcsp_instances, rcsp)
    return
end

