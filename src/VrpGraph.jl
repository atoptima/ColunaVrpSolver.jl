mutable struct VrpGraph{T}
    id::Int
    vert_ids::Vector{Int}
    bounds::Tuple{Float64,Float64}
    orig_sink::Int
    new_sink::Int
    mappings::Vector{Vector{VariableRef}}
    is_preproc::Bool
    model::T
    arcs::Vector{Tuple{Int,Int}}
    elem_sets::Vector{Vector{Int}}
    nb_resources::Int
    res_bounds::Vector{Vector{Tuple{Float64,Float64}}}
    res_is_main::Vector{Bool}
    res_is_binary::Vector{Bool}
    res_is_disposable::Vector{Bool}
    resid_to_binresid::Vector{Int}
    res_cons::Vector{Vector{Float64}}
    dist_matrix::Vector{Vector{Float64}}
    ng_sets::Vector{Vector{Int}}    # ng_sets[vertex_id] contains the elementarity sets for vertex_id
    src_id::Int
    snk_id::Int
    arc_ids::Vector{Cint}
end

mutable struct RCSPProblem
    graph::VrpGraph
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
    model::T, vertices::Vector{Int}, source::Int, sink::Int, bounds::Tuple{Int,Int},
) where T
    # For RCSP
    new_sink = sink
    vertices_ = copy(vertices)
    if sink == source
        new_sink = maximum(vertices) + 1
        push!(vertices_, new_sink)
    end

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
        model.nb_subproblems + 1, vert_ids, Float64.(bounds), sink, new_sink, Vector{VariableRef}[],
        false, model, Tuple{Int,Int}[], Vector{Int}[], 0, [Tuple{Float64,Float64}[] for _ in vertices_],
        Bool[], Bool[], Bool[], Int[], Vector{Float64}[], Vector{Float64}[], [Int[] for _ in eachindex(vert_ids)],
        src_id, snk_id, Cint[],
    )
    model.nb_subproblems += 1

    # cache the model's objective function for performance
    if model.form_obj == 0
        model.form_obj = objective_function(model.formulation)
    end
    return graph
end

function add_resource!(graph::VrpGraph; main=false, binary=false, disposable=true)
    if binary && disposable
        @error "Disposable binary resources are not supported"
    end
    id = graph.nb_resources
    graph.nb_resources += 1
    for b in graph.res_bounds
        push!(b, (0.0, 0.0))
    end
    push!(graph.res_is_main, main)
    push!(graph.res_is_binary, binary)
    push!(graph.res_is_disposable, disposable)
    push!(graph.resid_to_binresid, ifelse(binary, count(graph.res_is_binary), -1))
    return id
end

set_resource_bounds!(graph::VrpGraph, vertid::Int, resid::Int, lb::Int, ub::Int) =
    set_resource_bounds!(graph, vertid, resid, Float64(lb), Float64(ub))
function set_resource_bounds!(
    graph::VrpGraph, vertid::Int, resid::Int, lb::Float64, ub::Float64,
)
    graph.res_bounds[graph.vert_ids[vertid+1]+1][resid+1] = (lb, ub)
    return
end

function add_arc!(graph::VrpGraph, tail::Int, head::Int)
    h = (head == graph.orig_sink) ? graph.new_sink : head
    id = length(graph.arcs)
    if id + 1 > length(graph.arcs)
        resize!(graph.arcs, id + 1)
    end
    graph.arcs[id+1] = (tail, h)
    push!(graph.res_cons, zeros(Float64, graph.nb_resources))    # used only for BaPCod
    return id
end

set_arc_consumption!(graph::VrpGraph, arcid::Int, resid::Int, cons::Int) =
    set_arc_consumption!(graph, arcid, resid, Float64(cons))
function set_arc_consumption!(graph::VrpGraph, arcid::Int, resid::Int, cons::Float64)
    graph.res_cons[arcid+1][resid+1] = cons
    return
end

function add_arc_var_mapping!(graph::VrpGraph{T}, arcid::Int, var::Vector{VariableRef}) where {T}
    for v in var
        add_arc_var_mapping!(graph, arcid, v)
    end
    return
end

function add_arc_var_mapping!(graph::VrpGraph{T}, arcid::Int, var::VariableRef) where {T}
    mapped = get_mappedvarids!(graph, arcid)
    push!(mapped, var)
    spids = get(graph.model.spids_by_var, var, zeros(Bool, graph.id))
    if haskey(graph.model.spids_by_var, var)
        old_length = length(spids)
        if graph.id > old_length
            resize!(spids, graph.id)
        end
        for i in (old_length+1):(graph.id-1)
            spids[i] = false
        end
        spids[graph.id] = true
    else
        spids[graph.id] = true
        graph.model.spids_by_var[var] = spids
    end
    return
end

function add_graph!(model::T, graph::VrpGraph) where {T<:AbstractVrpModel}
    if graph.id != length(model.rcsp_instances) + 1
        @error "Graphs should be added in order"
    end

    # Add the graph to the VRP model
    push!(model.rcsp_instances, RCSPProblem(graph))
    return
end

function set_packing_sets!(
    is_vertex::Bool, model::T, psets::Vector{Vector{Tuple{VrpGraph{T},Int}}},
) where {T<:AbstractVrpModel}
    model.is_vertex_psets = is_vertex
    model.packing_sets = [[(graph.id - 1, elemid) for (graph, elemid) in pset] for pset in psets]
    empty!(model.pset_to_id)
    for pset in psets
        first = fill(true, length(model.rcsp_instances))
        for (graph, elemid) in pset
            if first[graph.id]
                push!(graph.elem_sets, [elemid])
                first[graph.id] = false
            else
                push!(graph.elem_sets[end], elemid)
            end
        end
    end
    for pset in model.packing_sets
        model.pset_to_id[pset] = length(model.pset_to_id)
    end
end

function set_vertex_packing_sets!(
    model::T, psets::Vector{Vector{Tuple{VrpGraph{T},Int}}},
) where {T<:AbstractVrpModel}
    set_packing_sets!(true, model, psets)
end

function set_arc_packing_sets!(
    model::T, psets::Vector{Vector{Tuple{VrpGraph{T},Int}}},
) where {T<:AbstractVrpModel}
    set_packing_sets!(false, model, psets)
end

function define_elementarity_sets_distance_matrix!(
    ::T, graph::VrpGraph, distmatrix::Vector{Vector{Float64}},
) where {T<:AbstractVrpModel}
    # check if the vector of vectors distmatrix is a square matrix
    lengths = length.(distmatrix)
    nb_psets = length(lengths)
    if !all(lengths .== nb_psets)
        @error "Distance matrix is not a square matrix"
    end

    # set the distance matrix
    graph.dist_matrix = distmatrix
end

function add_elem_set_to_vertex_init_ng_neighbourhood!(
    ::T, graph::VrpGraph, vertex_id::Int, es_id::Int
) where {T<:AbstractVrpModel}
    if vertex_id < 0 || vertex_id >= length(graph.vert_ids) || graph.vert_ids[vertex_id+1] == -1
        @error "Unknown vertex $vertex_id"
    end
    if es_id <= 0 || es_id > length(graph.elem_sets)
        @error "Unknown elementarity set $es_id"
    end
    if !(es_id in graph.ng_sets[vertex_id+1])
        push!(graph.ng_sets[vertex_id+1], es_id)
    end
    return
end

function add_capacity_cut_separator!(
    model::M, demandsets::Vector{Tuple{Vector{Tuple{VrpGraph{M},Int}},Float64}},
    capacity::Float64,
) where {M<:AbstractVrpModel}
    # check that all demand sets are packing sets and map all graph vertices to them
    vid_to_pset = [[-1 for _ in 1:length(rcsp.graph.vert_ids)] for rcsp in model.rcsp_instances]
    dem_sets = [([(ps[1].id - 1, ps[1].vert_ids[ps[2]+1]) for ps in ps_set], d) for (ps_set, d) in demandsets]
    for (ps_set, _) in dem_sets
        psid = get(model.pset_to_id, ps_set, -1)
        (psid == -1) && error(
            "Collection that is not a packing set was used in a capacity cut separator." *
            " Only the packing set collections can be used for add_capacity_cut_separator",
        )
        for ps in ps_set
            vid_to_pset[ps[1]+1][ps[2]] = psid
        end
    end

    # create and map variables to all uncovered arcs connecting packing set pairs
    nb_psets = length(dem_sets)
    arcs_by_pset_pair = [Tuple{Int,Int}[] for _ in 1:nb_psets, _ in 1:nb_psets]
    for gid in eachindex(model.rcsp_instances)
        graph = model.rcsp_instances[gid].graph
        for (id, (h, t)) in enumerate(graph.arcs)
            if isempty(graph.mappings[id])
                head = graph.vert_ids[h+1] + 1
                tail = graph.vert_ids[t+1] + 1
                edge = (head < tail) ? (head, tail) : (tail, head)
                if (vid_to_pset[gid][head] != -1) && (vid_to_pset[gid][tail] != -1)
                    push!(
                        arcs_by_pset_pair[vid_to_pset[gid][edge[1]]+1, vid_to_pset[gid][edge[2]]+1],
                        (gid, id - 1),
                    )
                end
            end
        end
    end
    id_demands = [Cint(0) for _ in 1:length(model.packing_sets)]
    for (ps_set, d) in dem_sets
        ps_id = model.pset_to_id[ps_set]
        id_demands[ps_id+1] = Cint(d)
    end
    num_missing_arcs = 0
    uncovered = Tuple{Int,Int}[]
    dims_psp = size(arcs_by_pset_pair)
    for head in 1:dims_psp[1], tail in (head+1):dims_psp[2]
        if (id_demands[head] > 0) && (id_demands[tail] > 0) && !isempty(arcs_by_pset_pair[head, tail])
            push!(uncovered, (head, tail))
            num_missing_arcs += length(arcs_by_pset_pair[head, tail])
        end
    end
    if length(uncovered) > 0
        println("VrpSolver: adding $(length(uncovered)) internal variables mapping to ",
            "$num_missing_arcs arcs for use by capacity cuts",
        )
    end
    if num_missing_arcs > 0
        @variable(model.formulation,
            RCCsepX[ps_pair in uncovered], Int
        )
        for (head, tail) in uncovered
            for (gid, arcid) in arcs_by_pset_pair[head, tail]
                graph = model.rcsp_instances[gid].graph
                add_arc_var_mapping!(graph, arcid, RCCsepX[(head, tail)])
            end
        end
    end
    push!(model.rcc_demands, (Cint(capacity), id_demands))
end

