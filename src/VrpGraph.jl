mutable struct VrpGraph{T}
    id::Int
    cptr::Ptr{Cvoid}
    vert_ids::Vector{Int}
    bounds::Tuple{Float64, Float64}
    orig_sink::Int
    new_sink::Int
    mappings::Vector{Vector{VariableRef}}
    is_preproc::Bool
    model::T
    arcs::Vector{Tuple{Int, Int}}
    elem_sets::Vector{Vector{Int}}
    nb_resources::Int
    res_bounds::Vector{Vector{Tuple{Float64, Float64}}}
    res_cons::Vector{Vector{Float64}}
    dist_matrix::Vector{Vector{Float64}}
    src_id::Int
    snk_id::Int
    arc_ids::Vector{Cint}
end

function get_mappedvarids(g::VrpGraph, arcid::Int)
    if !isassigned(g.mappings, arcid + 1)
        return VariableRef[]
    end
    return g.mappings[arcid+1]
end

function get_mappedvarids!(g::VrpGraph, arcid::Int)
    if !isassigned(g.mappings, arcid + 1)
        if arcid >= length(g.mappings)
            resize!(g.mappings, arcid + 1)
        end
        g.mappings[arcid+1] = VariableRef[]
    end
    return g.mappings[arcid+1]
end

function VrpGraph(
    model::T, vertices::Vector{Int}, source::Int, sink::Int, bounds::Tuple{Int, Int},
) where T
    # For RCSP
    new_sink = sink
    vertices_ = copy(vertices)
    if sink == source
        new_sink = maximum(vertices) + 1
        push!(vertices_, new_sink)
    end
    cptr_ = @try_ccall(
        (:createGraph_c, rcsp_path), Ptr{Cvoid}, (Cint, Ptr{Cint}, Cint, Cint),
        Cint(length(vertices_)), [Cint(v) for v in vertices_], Cint(source), Cint(new_sink),
    )

    # For BaPCod
    vert_ids = fill(-1, maximum(vertices_) + 1)
    src_id = 0
    snk_id = 0
    for (i, v) in enumerate(vertices_)
        vert_ids[v+1] = i - 1
        if v == source
            src_id = i - 1
        end
        if v == new_sink
            snk_id = i - 1
        end
    end

    # Create the graph object
    graph = VrpGraph(
        model.nb_subproblems + 1, cptr_, vert_ids, Float64.(bounds), sink, new_sink, Vector{VariableRef}[],
        false, model, Tuple{Int, Int}[], Vector{Int}[], 0, [Tuple{Float64, Float64}[] for _ in vertices_],
        Vector{Float64}[], Vector{Float64}[], src_id, snk_id, Cint[],
    )
    model.nb_subproblems += 1

    # cache the model's objective function for performance
    if model.form_obj == 0
        model.form_obj = objective_function(model.formulation)
    end
    return graph
end

function add_resource!(graph::VrpGraph; main = false)
    if rcsp_path != ""
        id = Int(ccall((:addResource_c, rcsp_path), Cint, (Ptr{Cvoid}, Cint), graph.cptr, Cint(main)))
    else
        id = graph.nb_resources
    end
    graph.nb_resources += 1
    for b in graph.res_bounds
        push!(b, (0.0, 0.0))
    end
    return id
end

function set_resource_bounds!(
    graph::VrpGraph, vertid::Int, resid::Int, lb::Float64, ub::Float64,
)
    vid = (vertid == graph.orig_sink) ? graph.new_sink : vertid
    @try_ccall(
        (:setResourceBounds_c, rcsp_path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64, Float64),
        graph.cptr, Cint(vid), Cint(resid), lb, ub,
    )
    graph.res_bounds[graph.vert_ids[vertid+1]+1][resid+1] = (lb, ub)
    return
end

function add_arc!(graph::VrpGraph, tail::Int, head::Int)
    h = (head == graph.orig_sink) ? graph.new_sink : head
    if rcsp_path != ""
        id = Int(ccall(
            (:addArc_c, rcsp_path), Cint, (Ptr{Cvoid}, Cint, Cint), graph.cptr, Cint(tail), Cint(h),
        ))
    else
        id = length(graph.arcs)
    end
    if id + 1 > length(graph.arcs)
        resize!(graph.arcs, id + 1)
    end
    graph.arcs[id+1] = (tail, h)
    push!(graph.res_cons, zeros(Float64, graph.nb_resources))    # used only for BaPCod
    return id
end

function set_arc_consumption!(graph::VrpGraph, arcid::Int, resid::Int, cons::Float64)
    @try_ccall(
        (:setArcConsumption_c, rcsp_path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64),
        graph.cptr, Cint(arcid), Cint(resid), cons,
    )
    graph.res_cons[arcid+1][resid+1] = cons
    return
end

function add_arc_var_mapping!(graph::VrpGraph{T}, arcid::Int, var::VariableRef) where {T}
    varid = getvarid!(graph.model, var)
    cost = coefficient(graph.model.form_obj, var)
    @try_ccall(
        (:addArcVarMapping_c, rcsp_path), Cvoid, (Ptr{Cvoid}, Cint, Cint, Float64),
        graph.cptr, Cint(arcid), Cint(varid - 1), cost,
    )
    mapped = get_mappedvarids!(graph, arcid)
    push!(mapped, var)
    if haskey(graph.model.spids_by_var, var)
        push!(graph.model.spids_by_var[var], graph.id)
    else
        graph.model.spids_by_var[var] = [graph.id]
    end
    return
end

function add_graph!(model::T, graph::VrpGraph) where {T <: AbstractVrpModel}
    if graph.id != length(model.rcsp_instances) + 1
        @error "Graphs should be added in order"
    end

    # Add the graph to the VRP model
    push!(model.rcsp_instances, RCSPProblem(graph))
    return
end

function preprocess_graph!(graph::VrpGraph)
    if !graph.is_preproc
        @try_ccall((:preprocessGraph_c, rcsp_path), Cint, (Ptr{Cvoid},), graph.cptr)
        graph.is_preproc = true
    end
    return
end

function set_vertex_packing_sets!(
    model::T, psets::Vector{Vector{Tuple{VrpGraph{T}, Int}}},
) where {T <: AbstractVrpModel}
    sizes = Cint.(length.(psets))
    graphs = vcat([getfield.(getindex.(ps, 1), :cptr) for ps in psets]...)
    vertids = vcat([Cint.(getindex.(ps, 2)) for ps in psets]...)
    @try_ccall(
        (:setVertexPackingSets_c, rcsp_path), Cvoid,
        (Cint, Ptr{Cint}, Ref{Ptr{Cvoid}}, Ptr{Cint}),
        Cint(length(psets)), sizes, graphs, vertids,
    )
    model.packing_sets = [[(graph.id - 1, vertid) for (graph, vertid) in pset] for pset in psets]
    for pset in psets
        first = fill(true, length(model.rcsp_instances))
        for (graph, vertid) in pset
            if first[graph.id]
                push!(graph.elem_sets, [vertid])
                first[graph.id] = false
            else
                push!(graph.elem_sets[end], vertid)
            end
        end
    end
end

function define_elementarity_sets_distance_matrix!(
    ::T, graph::VrpGraph, distmatrix::Vector{Vector{Float64}},
) where {T <: AbstractVrpModel}
    # check if the vector of vectors distmatrix is a square matrix
    lengths = length.(distmatrix)
    nb_psets = length(lengths)
    if !all(lengths .== nb_psets)
        @error "Distance matrix is not a square matrix"
    end

    # set the distance matrix
    dists = vcat(distmatrix...)
    @try_ccall(
        (:defineElemSetsDistMatrix_c, rcsp_path), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Float64}),
        graph.cptr, nb_psets, dists,
    )
    graph.dist_matrix = distmatrix
end
