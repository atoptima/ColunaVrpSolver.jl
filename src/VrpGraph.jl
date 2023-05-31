# Get the binary library path
global path = "$(ENV["RCSP_LIB_PATH"])"

abstract type AbstractVrpModel end

mutable struct VrpGraph{T}
    id::Int
    cptr::Ptr{Cvoid}
    bounds::Tuple{Float64, Float64}
    orig_sink::Int
    new_sink::Int
    mappings::Vector{Vector{VariableRef}}
    is_preproc::Bool
    model::T
    arcs::Dict{Int,Tuple{Int, Int}}
end

function get_mappedvarids(g::VrpGraph, arcid::Int)
    if !isassigned(g.mappings, arcid + 1)
        return VariableRef[]
    end
    return g.mappings[arcid + 1]
end

function get_mappedvarids!(g::VrpGraph, arcid::Int)
    if !isassigned(g.mappings, arcid + 1)
        if arcid >= length(g.mappings)
            resize!(g.mappings, arcid + 1)
        end
        g.mappings[arcid + 1] = VariableRef[]
    end
    return g.mappings[arcid + 1]
end

function VrpGraph(
    model::T, vertices::Vector{Int}, source::Int, sink::Int, bounds::Tuple{Int, Int}
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
        0, cptr_, Float64.(bounds), sink, new_sink, Vector{VariableRef}[], false, model,
        Dict{Int,Tuple{Int, Int}}()
    )

    # cache the model's objective function for performance
    if model.form_obj == 0
        model.form_obj = objective_function(model.formulation)
    end
    return graph
end

function add_resource!(graph::VrpGraph; main = false)
    id = ccall((:addResource_c, path), Cint, (Ptr{Cvoid}, Cint), graph.cptr, Cint(main))
    return Int(id)
end

function set_resource_bounds!(
    graph::VrpGraph, vertid::Int, resid::Int, lb::Float64, ub::Float64
)
    vid = (vertid == graph.orig_sink) ? graph.new_sink : vertid
    ccall(
        (:setResourceBounds_c, path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64, Float64),
        graph.cptr, Cint(vid), Cint(resid), lb, ub
    )
    return
end

function add_arc!(graph::VrpGraph, tail::Int, head::Int)
    h = (head == graph.orig_sink) ? graph.new_sink : head
    id = Int(ccall(
        (:addArc_c, path), Cint, (Ptr{Cvoid}, Cint, Cint), graph.cptr, Cint(tail), Cint(h)
    ))
    graph.arcs[id] = (tail, head)
    return id
end

function set_arc_consumption!(graph::VrpGraph, arcid::Int, resid::Int, cons::Float64)
    ccall(
        (:setArcConsumption_c, path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64),
        graph.cptr, Cint(arcid), Cint(resid), cons
    )
    return
end

function add_arc_var_mapping!(graph::VrpGraph{T}, arcid::Int, var::VariableRef) where {T}
    varid = getvarid!(graph.model, var)
    cost = coefficient(graph.model.form_obj, var)
    ccall(
        (:addArcVarMapping_c, path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64),
        graph.cptr, Cint(arcid), Cint(varid - 1), cost
    )
    mapped = get_mappedvarids!(graph, arcid)
    push!(mapped, var)
    return
end

function add_graph!(model::T, graph::VrpGraph) where {T <: AbstractVrpModel}
    # Add the graph to the VRP model
    push!(model.rcsp_instances, RCSPProblem(graph))

    # Save the graph id and return
    graph.id = length(model.rcsp_instances)
    return
end

function preprocess_graph!(graph::VrpGraph)
    if !graph.is_preproc
        ccall((:preprocessGraph_c, path), Cint, (Ptr{Cvoid},), graph.cptr)
        graph.is_preproc = true
    end
    return
end

function set_vertex_packing_sets!(
    ::T, psets::Vector{Vector{Tuple{VrpGraph{T}, Int}}}
)  where {T <: AbstractVrpModel}
    sizes = Cint.(length.(psets))
    graphs = vcat([getfield.(getindex.(ps, 1), :cptr) for ps in psets]...)
    vertids = vcat([Cint.(getindex.(ps, 2)) for ps in psets]...)
    ccall(
        (:setVertexPackingSets_c, path), Cvoid,
        (Cint, Ptr{Cint}, Ref{Ptr{Cvoid}}, Ptr{Cint}),
        Cint(length(psets)), sizes, graphs, vertids
    )
end

function define_elementarity_sets_distance_matrix!(
    ::T, graph::VrpGraph, distmatrix::Vector{Vector{Float64}}
) where {T <: AbstractVrpModel}
    # check if the vector of vectors distmatrix is a square matrix
    lengths = length.(distmatrix)
    nb_psets = length(lengths)
    if !all(lengths .== nb_psets)
        @error "Distance matrix is not a square matrix"
    end

    # set the distance matrix
    dists = vcat(distmatrix...)
    ccall(
        (:defineElemSetsDistMatrix_c, path), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Float64}),
        graph.cptr, nb_psets, dists
    )
end
