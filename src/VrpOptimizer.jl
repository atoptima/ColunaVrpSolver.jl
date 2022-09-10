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
function solve_RCSP_pricing(cbdata::CB, model::VrpModel, rcosts::Vector{Float64}) where {CB}
    # Get the reduced costs
    @show model.variables_by_id
    for vid in 1:get_maxvarid(model)
        @show vid
        rcosts[vid] = BlockDecomposition.callback_reduced_cost(
            cbdata, model.variables_by_id[vid]
        )
    end
    # @show rcosts

    # call the pricing solver
    rcsp = model.rcsp_instances[BlockDecomposition.callback_spid(cbdata, model.formulation)]
    paths = run_rcsp_pricing(rcsp, 0, rcosts)
    # @show paths

    # submit the priced paths to Coluna
    min_rc = Inf
    for p in paths
        solvals = Float64[]
        solvars = VariableRef[]
        varcount = Dict{Int, Int}()
        rc = 0.0
        for a in p
            for v in rcsp.graph.mappings[a + 1]
                vid = model.varids_by_var[v]
                varcount[vid] = get(varcount, vid, 0) + 1
                rc += rcosts[vid]
            end
        end
        min_rc = min(rc, min_rc)
        for vid in keys(varcount)
            push!(solvals, Float64(varcount[vid]))
            push!(solvars, model.variables_by_id[vid])
        end
        MathOptInterface.submit(
            model.formulation, BlockDecomposition.PricingSolution(cbdata),
            rc, solvars, solvals, PathVarData(g.indice, p)
        )
    end
    MathOptInterface.submit(
        model.formulation, BlockDecomposition.PricingDualBound(cbdata),
        min_rc
    )
end

# Define the function to separate capacity cuts via RCSP library
function separate_capacity_cuts(cbdata::CBD, model::VrpModel) where {CBD}
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
        # for m in cut.members
        #     val = callback_value(cbdata, model.variables_by_id[m.varid])
        #     if val > 1e-5
        #         @show m.coeff, m.varid, val
        #     end
        # end
        # lhs = sum(
        #     m.coeff * callback_value(cbdata, model.variables_by_id[m.varid])
        #     for m in cut.members
        # )
        # @show lhs, domain(cut)
        moi_cut = ScalarConstraint(
            sum(
                m.coeff *  model.variables_by_id[m.varid] for m in cut.members
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

    # apply the decomposition and store the axis
    nb_graphs = length(model.rcsp_instances)
    @axis(VrpGraphs, 1:nb_graphs)
    @dantzig_wolfe_decomposition(model.formulation, decomp, VrpGraphs)
    model.bd_graphs[1] = decomp

    # Register the custom variable data
    customvars!(model.formulation, PathVarData)

    # Set the mapped variables as representatives
    subproblems = getsubproblems(decomp)
    for var in model.variables_by_id
        subproblemrepresentative(var, subproblems)
    end

    # set the capacity cut callback
    if !isempty(model.rcc_separators)
        MathOptInterface.set(
            model.formulation, MathOptInterface.UserCutCallback(),
            (cbdata -> separate_capacity_cuts(cbdata, model))
        )
    end

    # preallocate a vector to store the reduced costs
    rcosts = zeros(Float64, get_maxvarid(model))

    # set the solution multiplicities and the pricing callback for each graph
    for g in 1:nb_graphs
        (L, U) = model.rcsp_instances[g].graph.bounds
        specify!(
            subproblems[g], lower_multiplicity = L, upper_multiplicity = U,
            solver = (
                cbdata -> solve_RCSP_pricing(cbdata, model, rcosts)
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
