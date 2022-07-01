struct CapacityCutMember
    graphid::Int
    arcid::Int
    coeff::Float64
end

struct CapacityCut
    rhs::Float64
    sense::Symbol
    members::Vector{CapacityCutMember}
end

domain(cut::CapacityCut) = (
    (cut.sense == :<=) ?
        MathOptInterface.LessThan(cut.rhs) : MathOptInterface.GreaterThan(cut.rhs)
)

function add_capacity_cut_separator!(
    model::VrpModel, demandsets::Vector{Tuple{Vector{Tuple{VrpGraph,Int}}, Float64}},
    capacity::Float64
)
    # make sure that all graphs are preprocessed
    for rcsp in model.rcsp_instances
        preprocess_graph!(rcsp.graph)
    end
    
    # get only the first pair (G, i) of each packing set because they all should share the same
    # packing set id
    graphs = [ps[1][1].cptr for (ps, _) in demandsets]
    vertids = [Cint(ps[1][2]) for (ps, _) in demandsets]

    # call the C++ function to create the capacity cut separator and add it to the model
    demands = [d for (_, d) in demandsets]
    push!(
        model.rcc_separators,
        ccall(
            (:addCapacityCutSeparator_c, path), Ptr{Cvoid},
            (Cint, Ref{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Float64}, Float64),
            Cint(length(demandsets)), graphs, vertids, demands, capacity
        )
    )
end

function run_capacity_cut_separators(model::VrpModel, sol::VrpSolution)
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
            cut_arcids = Vector{Cint}(undef, bufsize[1])
            cut_coeffs = Vector{Float64}(undef, bufsize[1])
            ccall(
                (:getCapacityCuts_c, path), Cvoid,
                (
                    Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                    Ptr{Float64}
                ),
                sep, senses, rhss, cut_starts, cut_graphids, cut_arcids, cut_coeffs
            )

            # convert and store the cuts
            for c in 1:nb_cuts
                cut = CapacityCut(
                    rhss[c], (senses[c] == 0) ? :>= : :<=, CapacityCutMember[]
                )
                nextcut!(model.coeffmanager, model)
                for m in cut_starts[c]:(cut_starts[c+1] - 1)
                    g = cut_graphids[m + 1] + 1
                    a = Int(cut_arcids[m + 1])
                    if !hascoeff(model.coeffmanager, g, a)
                        push!(cut.members, CapacityCutMember(g, a, cut_coeffs[m + 1]))
                        regcoeff!(model.coeffmanager, g, a)
                    end
                end
                isempty(cut.members) && continue
                push!(cuts, cut)
            end
        end
    end
    return cuts
end
