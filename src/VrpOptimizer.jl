mutable struct VrpOptimizer
    model::VrpModel
    status_dict::Dict{MathOptInterface.TerminationStatusCode, Symbol}
end

struct PathVarData <: BlockDecomposition.AbstractCustomData
    graphid::Int
    arcids::Vector{Int}
end

struct VrpSolution
    paths::Vector{Tuple{Float64, PathVarData}}
end

# Define the function to perform pricing via RCSP library
function solve_RCSP_pricing(
    cbdata::CBDataType, model::VrpModel, A::Vector{Vector{Int}}, __arc::VarRefArrayType,
    arc_rcosts::Vector{Vector{Float64}}, ismapped::Vector{Vector{Bool}}
) where {CBDataType, VarRefArrayType}

    # Get the reduced costs of arc variables
    g = BlockDecomposition.callback_spid(cbdata, model.formulation)
    for a in A[g]
        arc_rcosts[g][a + 1] = BlockDecomposition.callback_reduced_cost(
            cbdata, __arc[g, a]
        )
    end
    setup_rc = Coluna.MathProg.getcurcost(cbdata.form, cbdata.form.duty_data.setup_var)
    # @show arc_rcosts

    # call the pricing solver
    paths = run_rcsp_pricing(model.rcsp_instances[g], 0, arc_rcosts[g])
    # @show paths

    # submit the priced paths to Coluna
    min_rc = Inf
    for p in paths
        solvals = Float64[]
        solvars = JuMP.VariableRef[]
        arccount = Dict{Int, Int}()
        rc = 0.0
        for a in p
            arccount[a] = get(arccount, a, 0) + 1
            rc += arc_rcosts[g][a + 1]
        end
        min_rc = min(rc, min_rc)
        for a in keys(arccount)
            if ismapped[g][a + 1]
                push!(solvals, Float64(arccount[a]))
                push!(solvars, __arc[g, a])
            end
        end
        MathOptInterface.submit(
            model.formulation, BlockDecomposition.PricingSolution(cbdata),
            rc, solvars, solvals, PathVarData(g.indice, p)
        )
    end
    MathOptInterface.submit(
        model.formulation, BlockDecomposition.PricingDualBound(cbdata),
        min_rc + setup_rc # FIXME: only when Coluna adds setup_rc internally
    )
end

# Define the function to separate capacity cuts via RCSP library
function separate_capacity_cuts(
    cbdata::CBDataType, model::VrpModel, __arc::VarRefArrayType
) where {CBDataType, VarRefArrayType}
    # get the solution of the current relaxation
    sol = VrpSolution(Tuple{Float64, PathVarData}[])
    for (varid, varval) in cbdata.orig_sol
        var = Coluna.MathProg.getvar(cbdata.form, varid)
        if typeof(var.custom_data) == PathVarData
            push!(sol.paths, (varval, var.custom_data))
        end
    end

    # call the separator provided by the library
    cuts = run_capacity_cut_separators(model, sol)

    # add the separated cuts
    for cut in cuts
        for m in cut.members
            val = callback_value(cbdata, __arc[m.graphid, m.arcid])
            # if val > 1e-5
            #     @show m.coeff, m.graphid, m.arcid, val
            # end
        end
        lhs = sum(
            m.coeff * callback_value(cbdata, __arc[m.graphid, m.arcid]) for m in cut.members
        )
        # @show lhs, domain(cut)
        moi_cut = ScalarConstraint(
            sum(
                m.coeff *  __arc[m.graphid, m.arcid] for m in cut.members
            ),
            domain(cut)
        )
        MathOptInterface.submit(
            model.formulation, MathOptInterface.UserCut(cbdata), moi_cut
        )
    end
end

function VrpOptimizer(model::VrpModel, _::String, _::AbstractString)
    build_solvers!(model)
    graphs = getfield.(model.rcsp_instances, :graph)

    # add the mapped arc-variables to the model
    A = [collect(keys(graphs[g].mappings)) for g in 1:length(graphs)]
    @axis(VrpGraphs, 1:length(graphs))
    @variable(model.formulation, __arc[g in VrpGraphs, a in A[g]] >= 0)

    # compute the inverse of the mapping function
    map_inverse = Dict{VariableRef, Vector{Tuple{Int, Int}}}()
    for graph in graphs
        for (arcid, mapped_modelvars) in graph.mappings
            for modelvar in mapped_modelvars
                mapped_arcs = get(map_inverse, modelvar, Tuple{Int, Int}[])
                if isempty(mapped_arcs)
                    map_inverse[modelvar] = mapped_arcs
                end
                push!(mapped_arcs, (graph.id, arcid))
            end
        end
    end

    # add the mapping constraints to the model
    # TODO when switching to implicit master:
    # - keep a Dict of VariableRefs to ids
    # - use "obj = objective_function(m, AffExpr)" and "coefficient(obj, x)" to get
    #   the variable costs (also used when they are mapped)
    # - inform the RCSP solver about mappings and variable costs
    # - Convert the arcs to mapped variables in the pricing callback
    # - Treat the arcids as varids in the capacity cut separation
    @constraint(
        model.formulation, __map[v in collect(keys(map_inverse))],
        v == sum(__arc[g, a] for (g, a) in map_inverse[v])
    )

    # apply the decomposition and store the axis
    @dantzig_wolfe_decomposition(model.formulation, decomp, VrpGraphs)
    model.bd_graphs[1] = decomp

    # Register the custom variable data
    BlockDecomposition.customvars!(model.formulation, PathVarData)

    # build a structure to quickly filter non-mapped arc ids
    ismapped = [zeros(Bool, graph.max_arcid + 1) for graph in graphs]
    for g in 1:length(A)
        for a in A[g]
            ismapped[g][a + 1] = true
        end
    end

    # set the capacity cut callback
    if !isempty(model.rcc_separators)
        MathOptInterface.set(
            model.formulation, MathOptInterface.UserCutCallback(),
            (cbdata -> separate_capacity_cuts(cbdata, model, __arc))
        )
    end

    # preallocate a vector to store the reduced costs
    arc_rcosts = [zeros(Float64, graph.max_arcid + 1) for graph in graphs]

    # set the solution multiplicities and the pricing callback for each graph
    subproblems = getsubproblems(decomp)
    for g in 1:length(graphs)
        (L, U) = graphs[g].bounds
        specify!(
            subproblems[g], lower_multiplicity = L, upper_multiplicity = U,
            solver = (
                cbdata -> solve_RCSP_pricing(cbdata, model, A, __arc, arc_rcosts, ismapped)
            )
        )
    end

    # create a dictionary of return values for compatibility with old applications
    sd = Dict(
        MathOptInterface.OPTIMAL => :Optimal,
        MathOptInterface.INFEASIBLE => :Infeasible,
        MathOptInterface.ITERATION_LIMIT => :UserLimit,
        MathOptInterface.TIME_LIMIT => :UserLimit,
        MathOptInterface.OPTIMIZE_NOT_CALLED => :NotSolved
    )

    # create the optimizer and return it
    return VrpOptimizer(model, sd)
end

function set_cutoff!(opt::VrpOptimizer, cutoffvalue::Float64)
    objectiveprimalbound!(opt.model.formulation, cutoffvalue)
end

function JuMP.optimize!(opt::VrpOptimizer)
    optimize!(opt.model.formulation)
    status = get(opt.status_dict, JuMP.termination_status(opt.model.formulation), :Error)
    hassol = (JuMP.result_count(opt.model.formulation) > 0)
    return status, hassol
end

function get_objective_value(opt::VrpOptimizer)
    return objective_value(opt.model.formulation)
end

function get_value(::VrpOptimizer, var::JuMP.VariableRef)
    return value(var)
end
