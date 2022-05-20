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

mutable struct RCSPProblem
    graph::VrpGraph
    solver::Ptr{Cvoid}
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
    # Add it to the VRP model
    push!(model.rcsp_instances, RCSPProblem(graph, Ptr{Cvoid}(0)))
    graph.id = length(model.rcsp_instances)
    return
end

function run_rcsp_pricing(rcsp::RCSPProblem, phase::Int, var_rcosts::Vector{Float64})
    # Call the RCSP pricing solver
    output = [Ptr{Cvoid}(0)]
    nb_arcs = [Cint(0)]
    nb_paths = Int(ccall(
        (:runPricing_c, path), Cint,
        (Ptr{Cvoid}, Cint, Cint, Ptr{Float64}, Ref{Ptr{Cvoid}}, Ref{Cint}),
        rcsp.solver, Cint(phase), Cint(length(var_rcosts)), var_rcosts,
        output, nb_arcs
    ))

    # get the output paths, each path as a vertor of arc ids
    starts = Vector{Cint}(undef, nb_paths + 1)
    arcs = Vector{Cint}(undef, Int(nb_arcs[1]))
    ccall(
        (:getOutputPaths_c, path), Cvoid, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}),
        output[1], starts, arcs
    )

    # convert and return them
    paths = Vector{Int}[]
    for p in 1:nb_paths
        push!(paths, [Int(a) for a in arcs[(starts[p] + 1):starts[p + 1]]])
    end
    return paths
end

function set_vertex_packing_sets!(
    model::T, psets::Vector{Vector{Tuple{VrpGraph, Int}}}
)  where {T <: AbstractVrpModel}
    sizes = Cint.(length.(psets))
    graphs = vcat([getfield.(getindex.(ps, 1), :cptr) for ps in psets]...)
    vertids = vcat([Cint.(getindex.(ps, 2)) for ps in psets]...)
    ccall(
        (:setVertexPackingSets_c, path), Cvoid,
        (Cint, Ptr{Cint}, Ref{Ptr{Cvoid}}, Ptr{Cint}),
        Cint(length(psets)), sizes, graphs, vertids
    )
    # TODO: save a dictionaty of (graph, vertex) -> psetId in model for CC separation
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

function add_capacity_cut_separator!(
    model::T, demandsets::Vector{Tuple{Vector{Tuple{VrpGraph,Int}}, Float64}}, capacity::Float64
) where {T <: AbstractVrpModel}
    # get only the first pair (G, i) of each packing set because they all should share the same
    # packing set id
    graphs = [ps[1][1].cptr for (ps, _) in demandsets]
    vertids = [Cint(ps[1][2]) for (ps, _) in demandsets]

    # call the C++ function to create the capacity cut separator and add it to the model
    demands = [d for (_, d) in demandsets]
    model.rcc_separator = ccall(
        (:addCapacityCutSeparator_c, path), Ptr{Cvoid},
        (Cint, Ref{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Float64}, Float64),
        Cint(length(demandsets)), graphs, vertids, demands, capacity
    )
end

function build_solvers!(model::T) where {T <: AbstractVrpModel}
    # Instantiate an RCSP solver for each RCSP problem instance
    for prob in model.rcsp_instances
        prob.solver = ccall(
            (:createAndPrepareSolver_c, path), Ptr{Cvoid}, (Ptr{Cvoid},),
            prob.graph.cptr
        )
    end
    return
end
