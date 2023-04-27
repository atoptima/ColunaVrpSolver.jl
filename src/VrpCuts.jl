abstract type AbstractVrpSolution end

struct CapacityCutMember
    varid::Int
    coeff::Float64
end

struct CapacityCut
    rhs::Float64
    sense::Symbol
    members::Vector{CapacityCutMember}
end

struct RankOneCutData <: BlockDecomposition.AbstractCustomData
    rhs::Float64
    data_ptr::Ptr{Cvoid}
    separator_ptr::Ptr{Cvoid}
end

struct RCCPreSeparator
    graphs::Vector{Ptr{Cvoid}}
    vertids::Vector{Cint}
    demands::Vector{Float64}
    capacity::Float64
end

domain(cut::CapacityCut) = (
    (cut.sense == :<=) ?
        MathOptInterface.LessThan(cut.rhs) : MathOptInterface.GreaterThan(cut.rhs)
)

function add_capacity_cut_separator!(
    model::M, demandsets::Vector{Tuple{Vector{Tuple{VrpGraph{M},Int}}, Float64}},
    capacity::Float64
) where {M <: AbstractVrpModel}
    # make sure that all graphs are preprocessed
    for rcsp in model.rcsp_instances
        preprocess_graph!(rcsp.graph)
    end
    
    # get only the first pair (G, i) of each packing set because they all should share the same
    # packing set id, which is the id used inside the RCSP C++ module to assign demands.
    # Note: Fixing this would require to change the VrpSolver user interface, making it less
    #       intuitive. In the current VrpSolver interface the same vertex grouping as the one
    #       defined by the packing sets is repeated in the argument demandsets to avoid working
    #       with packing set ids, and mimic the mathematical notation used in papers.
    graphs = [ps[1][1].cptr for (ps, _) in demandsets]
    vertids = [Cint(ps[1][2]) for (ps, _) in demandsets]

    # call the C++ function to create the capacity cut separator and add it to the model
    demands = [d for (_, d) in demandsets]
    push!(
        model.rcc_pre_separators, RCCPreSeparator(graphs, vertids, demands, capacity)
    )
end

function build_capacity_cut_separators!(model::M) where {M <: AbstractVrpModel}
    for presep in model.rcc_pre_separators
        push!(
            model.rcc_separators,
            ccall(
                (:addCapacityCutSeparator_c, path), Ptr{Cvoid},
                (Cint, Ref{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Float64}, Float64, Ptr{Cvoid}),
                Cint(length(presep.demands)), presep.graphs, presep.vertids, presep.demands,
                presep.capacity,
                get_rcsp_params(model.parameters[1], PARAM_CLASS_ROUND_CAP_CUTS_SEPARATOR)
            )
        )
    end
end

function run_capacity_cut_separators(
    model::M, sol::S
) where {M <: AbstractVrpModel, S <: AbstractVrpSolution}
    # prepare the fractional solution data
    path_values = Float64[]
    path_graphids = Cint[]
    path_starts = Cint[]
    arcids = Cint[]
    bufpos = 1
    for (val, path) in sol.paths
        push!(path_values, val)
        push!(path_graphids, path.graphid - 1)
        push!(path_starts, bufpos - 1)
        for a in path.arcids
            push!(arcids, a)
            bufpos += 1
        end
    end
    push!(path_starts, bufpos - 1)

    # accumulate the capacity cuts for all defined separators
    cuts = CapacityCut[]
    for sep in model.rcc_separators
        # call the separator and get the output size
        bufsize = [Cint(0)]
        nb_cuts = ccall(
            (:separateCapacityCuts_c, path), Cint,
            (Ptr{Cvoid}, Cint, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
            sep, Cint(length(sol.paths)), path_values, path_graphids, path_starts, arcids,
            bufsize
        )

        if nb_cuts > 0
            # allocate the buffer and retrieve the cuts
            senses = Vector{Cint}(undef, nb_cuts)
            rhss = Vector{Cint}(undef, nb_cuts)
            cut_starts = Vector{Cint}(undef, nb_cuts + 1)
            cut_graphids = Vector{Cint}(undef, bufsize[1])
            cut_varids = Vector{Cint}(undef, bufsize[1])
            cut_coeffs = Vector{Float64}(undef, bufsize[1])
            ccall(
                (:getCapacityCuts_c, path), Cvoid,
                (
                    Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                    Ptr{Float64}
                ),
                sep, senses, rhss, cut_starts, cut_graphids, cut_varids, cut_coeffs
            )

            # convert and store the cuts
            for c in 1:nb_cuts
                cut = CapacityCut(
                    rhss[c], (senses[c] == 0) ? :>= : :<=, CapacityCutMember[]
                )
                nextcut!(model.coeffmanager, model)
                for m in cut_starts[c]:(cut_starts[c+1] - 1)
                    # g = cut_graphids[m + 1] + 1
                    v = Int(cut_varids[m + 1] + 1)
                    if !hascoeff(model.coeffmanager, v)
                        push!(cut.members, CapacityCutMember(v, cut_coeffs[m + 1]))
                        regcoeff!(model.coeffmanager, v)
                    end
                end
                isempty(cut.members) && continue
                push!(cuts, cut)
            end
        end
    end
    return cuts
end

function build_rank_one_cut_separator!(model::M) where {M <: AbstractVrpModel}
    # make sure that all graphs are preprocessed
    for rcsp in model.rcsp_instances
        preprocess_graph!(rcsp.graph)
    end

    # call the C++ code to build the separator
    model.rank1cut_separator = ccall(
        (:addRankOneCutSeparator_c, path), Ptr{Cvoid}, (Cint, Ref{Ptr{Cvoid}}, Ptr{Cvoid}),
        Cint(length(model.rcsp_instances)), map(x -> x.graph.cptr, model.rcsp_instances),
        get_rcsp_params(model.parameters[1], PARAM_CLASS_LIM_MEM_RANK_ONE_CUTS_SEPARATOR)
    )
end

function get_cutbuf_size(cutsep_phase::Int, max_cuts::Int)
    bufsize = max_cuts
    if cutsep_phase >= 3
        bufsize += div(max_cuts * 3, 2)
    end
    if cutsep_phase >= 4
        bufsize += max_cuts * (cutsep_phase - 3)
    end
    # @show bufsize
    return bufsize
end

function get_rank1cut_rhs(cutptr::Ptr{Cvoid})
    # @show cutptr
    return ccall((:getRankOneCutRhs_c, path), Float64, (Ptr{Cvoid},), cutptr)
end

function run_rank_one_cut_separator(
    model::M, sol::S
) where {M <: AbstractVrpModel, S <: AbstractVrpSolution}
    # prepare the fractional solution data
    path_values = Float64[]
    path_graphids = Cint[]
    path_starts = Cint[]
    arcids = Cint[]
    bufpos = 1
    for (val, path) in sol.paths
        push!(path_values, val)
        push!(path_graphids, path.graphid - 1)
        push!(path_starts, bufpos - 1)
        for a in path.arcids
            push!(arcids, a)
            bufpos += 1
        end
    end
    push!(path_starts, bufpos - 1)

    # call the separator to get the violated cut pointers
    max_cuts = get_rcsp_rank1cut_param_value(Int, model.parameters[1], :maxNumPerRound)
    mem_type = get_rcsp_rank1cut_param_value(Int, model.parameters[1], :memoryType)
    phase = [Cint(model.cutsep_phase)]
    cutbuf = Vector{Ptr{Cvoid}}(undef, get_cutbuf_size(model.cutsep_phase, max_cuts))
    nb_rank1cuts = ccall(
        (:separateRankOneCutCuts_c, path), Cint,
        (
            Ptr{Cvoid}, Ptr{Cint}, Cint, Cint, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
            Ref{Ptr{Cvoid}}
        ),
        model.rank1cut_separator, phase, mem_type, Cint(length(sol.paths)), path_values,
        path_graphids, path_starts, arcids, cutbuf
    )
    @show nb_rank1cuts
    model.cutsep_phase = phase[1]
    resize!(cutbuf, nb_rank1cuts)

    # convert and return the cuts
    cuts = [
        RankOneCutData(get_rank1cut_rhs(cutptr), cutptr, model.rank1cut_separator)
        for cutptr in cutbuf
    ]
    # rhss = getfield.(cuts, :rhs)
    # @show rhss
    # viols = [compute_violation(c, sol) for c in cuts]
    # @show viols
    return cuts
end

# function to compute rank-1 cut coefficients called by Coluna
function Coluna.MathProg.computecoeff(
    ::Coluna.Variable, var_custom_data::PathVarData,
    ::Coluna.Constraint, constr_custom_data::RankOneCutData
)
    return compute_coeff_from_data(var_custom_data, constr_custom_data)
end
function compute_coeff_from_data(
    var_custom_data::PathVarData, constr_custom_data::RankOneCutData
)
    arcids = var_custom_data.arcids
    return ccall(
        (:getRankOneCutCoeff_c, path), Float64, (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cvoid}),
        constr_custom_data.separator_ptr, Cint(var_custom_data.graphid - 1), Cint(length(arcids)),
        arcids, constr_custom_data.data_ptr
    )
end

# Compute the violation of a rank-1 cut for a given solution
function compute_violation(cut::RankOneCutData, sol::S) where {S <: AbstractVrpSolution}
    viol = -cut.rhs
    for (val, path) in sol.paths
        coeff = compute_coeff_from_data(path, cut)
        viol += coeff * val
    end
    return viol
end
