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

function Coluna.Algorithm.run!(
    algo::RedCostFixAndEnumAlgorithm, ::Coluna.Env,
    reform::Coluna.MathProg.Reformulation, input::Coluna.Algorithm.OptimizationState
)
    masterform = Coluna.MathProg.getmaster(reform)
    changed = false
    for (_, spform) in Coluna.MathProg.get_dw_pricing_sps(reform)
        cbdata = Coluna.MathProg.PricingCallbackData(spform)
        changed |= algo.func(masterform, cbdata, input)
    end
    return changed
end

function Coluna.Algorithm.run!(
    algo::SolveByMipAlgorithm, ::Coluna.Env, reform::Coluna.MathProg.Reformulation,
    input::Coluna.Algorithm.OptimizationState
)
    # Call the solver by MIP function
    masterform = Coluna.MathProg.getmaster(reform)
    return algo.func(masterform, input)
end

# Define the function to perform reduced-cost fixing and enumeration via RCSP library
function run_redcostfixing_and_enumeration!(
    masterform::Coluna.MathProg.Formulation, cbdata::CB, model::VrpModel, rcosts::Vector{Float64},
    optstate::Coluna.Algorithm.OptimizationState
) where {CB}
    # Get the primal and dual bounds
    dual_bnd = Coluna.Algorithm.get_lp_dual_bound(optstate).value
    primal_bnd = Coluna.Algorithm.get_ip_primal_bound(optstate).value

    # Get the dual variables associated to the subproblem bounds constraints
    dualsol = Coluna.Algorithm.get_best_lp_dual_sol(optstate)
    convdual = 0.0
    spid = BlockDecomposition.callback_spid(cbdata, model.formulation)
    reform = masterform.parent_formulation
    constrid = Coluna.MathProg.get_dw_pricing_sp_ub_constrid(reform, spid)
    convdual -= dualsol[constrid]
    constrid = Coluna.MathProg.get_dw_pricing_sp_lb_constrid(reform, spid)
    convdual -= dualsol[constrid]

    # Call the reduced cost fixing and enumeration function
    rcsp = model.rcsp_instances[spid]
    run_rcsp_rcostfix_and_enum(rcsp, rcosts, primal_bnd - dual_bnd - convdual)

    return false
end

# Define the function to solve the problem by MIP to be called in the node finalizer
function run_solve_by_mip!(
    masterform::Coluna.MathProg.Formulation, model::VrpModel,
    optstate::Coluna.Algorithm.OptimizationState
)
    # TODO: solve the problem of restoring this recorded RCSP state only once in the
    #       RCSP pricing
    
    # Record the RCSP princing state and store it
    storage = Coluna.MathProg.getstorage(masterform)
    unit = storage.units[VrpNodeInfoUnit].storage_unit
    for rcsp in model.rcsp_instances
        record_rcsp_state(rcsp)
    end
    unit.rcsp_state = [rcsp.state for rcsp in model.rcsp_instances]

    # Return an unchanged optimization state
    return Coluna.Algorithm.OptimizationState(
        masterform, 
        ip_primal_bound = Coluna.Algorithm.get_ip_primal_bound(optstate),
        termination_status = Coluna.OTHER_LIMIT
    )
end

# Define the function to perform pricing via RCSP library
function solve_RCSP_pricing(
    cbdata::CB, stage::Int, model::VrpModel, rcosts::Vector{Float64}
) where {CB}
    # Get the reduced costs
    for vid in 1:get_maxvarid(model)
        rcosts[vid] = BlockDecomposition.callback_reduced_cost(
            cbdata, model.variables_by_id[vid]
        )
    end
    setup_var = Coluna.MathProg.getname(cbdata.form, cbdata.form.duty_data.setup_var)
    @show setup_var
    var_rcosts = tuple.(model.variables_by_id, rcosts)
    @show cbdata.form.parent_formulation.lp_count
    @show var_rcosts

    # call the pricing solver
    rcsp = model.rcsp_instances[BlockDecomposition.callback_spid(cbdata, model.formulation)]
    paths = run_rcsp_pricing(rcsp, stage - 1, rcosts)
    if !isempty(paths)
        print("<st=$(stage - 1)>")  # TODO: move this to the Coluna's log
    end
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
            rc, solvars, solvals, PathVarData(rcsp.graph.id, p)
        )
    end
    MathOptInterface.submit(
        model.formulation, BlockDecomposition.PricingDualBound(cbdata),
        (stage == 1) ? min_rc : -Inf
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
            solver = [
                cbdata -> solve_RCSP_pricing(cbdata, stage, model, rcosts) for stage in 1:3
            ]
        )
    end

    # set the callback for reduced-cost fixing and enumeration, and solver by MIP
    model.redcostfix_enum_algo.func =
        (masterform, cbdata, optstate) -> run_redcostfixing_and_enumeration!(
            masterform, cbdata, model, rcosts, optstate
        )
    model.solve_by_mip_algo.func =
        (masterform, optstate) -> run_solve_by_mip!(masterform, model, optstate)

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
