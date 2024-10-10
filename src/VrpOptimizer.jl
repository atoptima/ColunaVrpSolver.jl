mutable struct VrpOptimizer
    model::VrpModel
    status_dict::Dict{MathOptInterface.TerminationStatusCode, Symbol}
end

struct VrpSolution <: AbstractVrpSolution
    paths::Vector{Tuple{Float64, PathVarData}}
end

function Coluna.Algorithm.run!(
    algo::RedCostFixAndEnumAlgorithm, ::Coluna.Env,
    reform::Coluna.MathProg.Reformulation, input::Coluna.Algorithm.OptimizationState,
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
    input::Coluna.Algorithm.OptimizationState,
)
    # Call the solver by MIP function
    masterform = Coluna.MathProg.getmaster(reform)
    pricing_sps = Coluna.MathProg.get_dw_pricing_sps(reform)
    cbdata_vec = Coluna.MathProg.PricingCallbackData[]
    for (_, spform) in pricing_sps
        push!(cbdata_vec, Coluna.MathProg.PricingCallbackData(spform))
    end
    return algo.func(masterform, cbdata_vec, input)
end

function should_solve_by_mip(unit::VrpNodeInfoUnit, model::VrpModel)
    if all(unit.enumerated)
        nb_paths = sum(get_number_of_enum_paths(rcsp) for rcsp in model.rcsp_instances)
        return (nb_paths <= model.parameters[1].coluna_vrp_params.RCSPmaxNumOfEnumSolutionsForMIP)
    end
    return false
end

function update_cutsep_status!(
    unit::VrpNodeInfoUnit, model::VrpModel, dual_bnd::Float64, curr_gap::Float64,
)
    # compute the gap ratio and update the last cut-round gap
    gap_diff_ratio = (unit.last_cutrnd_gap - curr_gap) / unit.last_cutrnd_gap
    if unit.separated_cuts
        println("Gap improvement since the last cut separation : $(gap_diff_ratio) ($(dual_bnd)))")
    end
    threshold = get_parameter(model, :CutTailingOffThreshold)
    # @show unit.last_cutrnd_gap, curr_gap, gap_diff_ratio, threshold
    unit.last_cutrnd_gap = curr_gap
    # @show unit.cutsep_phase
    if unit.cutsep_phase >= 0
        model.cutsep_phase = unit.cutsep_phase
        unit.cutsep_phase = -1
    end

    # Increase the tailing off counter if the gap did not reduce enough
    if unit.separated_cuts && gap_diff_ratio < threshold
        counter_threshold = get_parameter(model, :CutTailingOffCounterThreshold)
        unit.tailoff_counter += 1
        if unit.tailoff_counter < counter_threshold
            println(
                "Cut generation tailing off counter increased to $(unit.tailoff_counter)",
            )
            max_rows = get_rcsp_rank1cut_param_value(Int, model.parameters[1], :maxNumRows)
            if unit.tailoff_counter == (counter_threshold - 1) && model.cutsep_phase == 0
                model.cutsep_phase = min(3, max_rows)
            end
        else
            # set the flag to stop cut generation by tailing off
            unit.should_stop_cutsep = true
        end
    end
    unit.separated_cuts = false
    return nothing
end

# function to get the dual values of all rank-one cuts and clean-up the inactive ones
function get_rankonecut_duals!(cbdata::CB, cleanup::Bool) where {CB}
    cutdata = Ptr{Cvoid}[]
    duals = Float64[]
    masterform = cbdata.form.parent_formulation
    changed_form = false
    for (cid, constr) in Coluna.MathProg.getconstrs(masterform)
        if Coluna.MathProg.iscuractive(masterform, cid) &&
           Coluna.MathProg.getduty(cid) <= Coluna.MathProg.MasterUserCutConstr
            dualval = Coluna.MathProg.getcurincval(masterform, constr)
            if abs(dualval) > 1e-7 # FIXME: make 1e-7 configurable
                if typeof(constr.custom_data) == RankOneCutData
                    # ctrname = Coluna.MathProg.getname(masterform, constr)
                    # @show ctrname, dualval, constr.custom_data.data_ptr
                    push!(cutdata, constr.custom_data.data_ptr)
                    push!(duals, dualval)
                end
            elseif cleanup
                Coluna.MathProg.deactivate!(masterform, cid)
                changed_form = true
            end
        end
    end
    return cutdata, duals, changed_form
end

# Define the function to perform reduced-cost fixing and enumeration via RCSP library
function run_redcostfixing_and_enumeration!(
    masterform::Coluna.MathProg.Formulation, cbdata::CB, model::VrpModel, rcosts::Vector{Float64},
    optstate::Coluna.Algorithm.OptimizationState,
) where {CB}
    # Get the rank-1 cuts and their dual values
    r1cut_ptrs, r1cut_duals, changed_form = get_rankonecut_duals!(cbdata, true)

    # @info "In run_redcostfixing_and_enumeration!"
    storage = Coluna.MathProg.getstorage(masterform)
    unit = storage.units[VrpNodeInfoUnit].storage_unit
    if should_solve_by_mip(unit, model)
        return false
    end

    # Get the reduced costs
    for vid in 1:get_maxvarid(model)
        rcosts[vid] = BlockDecomposition.callback_reduced_cost(
            cbdata, model.variables_by_id[vid],
        )
    end
    # @show rcosts

    # Get the primal and dual bounds and compute the current gap
    dual_bnd = Coluna.Algorithm.get_lp_dual_bound(optstate).value
    primal_bnd = Coluna.Algorithm.get_ip_primal_bound(optstate).value
    curr_gap = (primal_bnd - dual_bnd) * 100.0 / primal_bnd

    # Update the cut separation phase
    update_cutsep_status!(unit, model, dual_bnd, curr_gap)

    # Stop if the gap did not reduce enough
    gap_ratio = curr_gap / unit.last_rcost_fix_gap
    threshold = get_parameter(model, :ReducedCostFixingThreshold)
    # @show unit.last_rcost_fix_gap, curr_gap, gap_ratio, threshold
    if gap_ratio > threshold
        println("Full reduced cost fixing is not called (gap ratio is $(gap_ratio))")
        return false
    end

    # Get the dual variables associated to the subproblem bounds constraints
    dualsol = Coluna.Algorithm.get_best_lp_dual_sol(optstate)
    convdual = 0.0
    reform = masterform.parent_formulation
    sp_uid = Coluna.MathProg.getuid(cbdata.form)
    constrid = Coluna.MathProg.get_dw_pricing_sp_ub_constrid(reform, sp_uid)
    convdual -= dualsol[constrid]
    constrid = Coluna.MathProg.get_dw_pricing_sp_lb_constrid(reform, sp_uid)
    convdual -= dualsol[constrid]

    # Call the reduced cost fixing and enumeration function
    spid = BlockDecomposition.callback_spid(cbdata, model.formulation)
    rcsp = model.rcsp_instances[spid]
    unit.enumerated[spid] |= run_rcsp_rcostfix_and_enum(
        rcsp, rcosts, r1cut_ptrs, r1cut_duals, primal_bnd - dual_bnd - convdual,
    )
    unit.last_rcost_fix_gap = curr_gap

    # deactivate the irrelevant columns from the master if enumerated
    changed = false
    if unit.enumerated[spid]
        colpaths = PathVarData[]
        for (vid, var) in Coluna.MathProg.getvars(masterform)
            if Coluna.MathProg.iscuractive(masterform, vid) &&
               Coluna.MathProg.getduty(vid) <= Coluna.MathProg.MasterCol &&
               var.custom_data.graphid == rcsp.graph.id
                push!(colpaths, var.custom_data)
            end
        end
        # @show colpaths
        isrelevant = check_enumerated_paths(rcsp, colpaths)
        # @show isrelevant
        col = 1
        for (vid, var) in Coluna.MathProg.getvars(masterform)
            if Coluna.MathProg.iscuractive(masterform, vid) &&
               Coluna.MathProg.getduty(vid) <= Coluna.MathProg.MasterCol &&
               var.custom_data.graphid == rcsp.graph.id
                # varname = Coluna.MathProg.getname(masterform, var)
                if !isrelevant[col]
                    Coluna.MathProg.deactivate!(masterform, vid)
                    changed = true
                end
                col += 1
            end
        end
    end
    # FIXME: if changed, solve the master again

    return false
end

# Define the function to solve the problem by MIP to be called in the node finalizer
function run_solve_by_mip!(
    masterform::Coluna.MathProg.Formulation, model::VrpModel, cbdata_vec::Vector{CB},
    optstate::Coluna.Algorithm.OptimizationState,
) where {CB}
    # Sort `cbdata_vec` by subproblem id (FIXME: will not work if the user creates other subproblems)
    sort!(cbdata_vec, by = cbd -> BlockDecomposition.callback_spid(cbd, model.formulation))

    # Check if should solve the problem by MIP
    storage = Coluna.MathProg.getstorage(masterform)
    unit = storage.units[VrpNodeInfoUnit].storage_unit
    if should_solve_by_mip(unit, model)
        # Solve the problem by MIP
        primal_bnd = Coluna.Algorithm.get_ip_primal_bound(optstate).value
        @time sol = solve_vrp_by_mip(masterform, model, cbdata_vec, primal_bnd)

        # Set the optimization state to output
        output = Coluna.Algorithm.OptimizationState(
            masterform,
            ip_primal_bound = Coluna.Algorithm.get_ip_primal_bound(optstate),
            termination_status = Coluna.OPTIMAL,
        )

        # Build the primal solution if any
        if !isempty(sol)
            Coluna.Algorithm.add_ip_primal_sol!(output, sol[1])
        end

        # Return the optimization state
        return output
    end

    # Record the RCSP princing state and store it
    unit.rcsp_states = [
        RCSPState(rcsp.solver, record_rcsp_state(rcsp)) for rcsp in model.rcsp_instances
    ]
    # @info "In run_solve_by_mip! $(unit.rcsp_states)"

    # Save the cut separation phase
    # @show model.cutsep_phase
    unit.cutsep_phase = model.cutsep_phase

    # Return an unchanged optimization state
    return Coluna.Algorithm.OptimizationState(
        masterform,
        ip_primal_bound = Coluna.Algorithm.get_ip_primal_bound(optstate),
        termination_status = Coluna.OTHER_LIMIT,
    )
end

# Define the function to perform pricing via RCSP library
function solve_RCSP_pricing(
    cbdata::CB, stage::Int, model::VrpModel, rcosts::Vector{Float64},
) where {CB}
    # Get the reduced costs and the non-robust cuts' duals
    for vid in 1:get_maxvarid(model)
        rcosts[vid] = BlockDecomposition.callback_reduced_cost(
            cbdata, model.variables_by_id[vid],
        )
    end
    r1cut_ptrs, r1cut_duals, _ = get_rankonecut_duals!(cbdata, false)
    # @show r1cut_ptrs
    # @show r1cut_duals
    # setup_var = Coluna.MathProg.getname(cbdata.form, cbdata.form.duty_data.setup_var)
    # @show setup_var
    # var_rcosts = tuple.(model.variables_by_id, rcosts)
    # @show cbdata.form.parent_formulation.lp_count
    # @show var_rcosts

    # call the pricing solver
    spid = BlockDecomposition.callback_spid(cbdata, model.formulation)
    rcsp = model.rcsp_instances[spid]
    masterform = cbdata.form.parent_formulation
    storage = Coluna.MathProg.getstorage(masterform)
    unit = storage.units[VrpNodeInfoUnit].storage_unit
    if isempty(unit.enumerated)
        unit.enumerated = zeros(Bool, length(model.rcsp_instances))
    end
    # _stage = unit.enumerated[spid] ? 0 : (stage - 1) # FIXME: Coluna needs solutions!
    paths = run_rcsp_pricing(rcsp, stage - 1, rcosts, r1cut_ptrs, r1cut_duals)
    # if unit.enumerated[spid]
    #     @show stage, paths
    # end
    # @show paths

    # submit the priced paths to Coluna
    min_rc = Inf
    for p in paths
        solvals = Float64[]
        solvars = VariableRef[]
        varcount = Dict{Int, Int}()
        rc = 0.0
        for a in p
            for v in rcsp.graph.mappings[a+1]
                vid = model.varids_by_var[v]
                varcount[vid] = get(varcount, vid, 0) + 1
                rc += rcosts[vid]
            end
        end
        path_data = PathVarData(rcsp.graph.id, p)
        arcs = [model.rcsp_instances[1].graph.arcs[a+1] for a in p]
        degrees = [length([a for a in arcs if a[2] == i]) for i in (1, 3, 5, 8, 20)]
        for i in eachindex(r1cut_ptrs)
            coeff = compute_coeff_from_data(
                path_data, RankOneCutData(0.0, r1cut_ptrs[i], model.rank1cut_separator),
            )
            # if coeff != 0.0 && length(arcs) == 8 && degrees == [1, 2, 1, 2, 1] 
            #     @show coeff, r1cut_duals[i]
            # end
            rc -= coeff * r1cut_duals[i]
        end
        # @show rc, arcs
        min_rc = min(rc, min_rc)
        for vid in keys(varcount)
            push!(solvals, Float64(varcount[vid]))
            push!(solvars, model.variables_by_id[vid])
        end
        MathOptInterface.submit(
            model.formulation, BlockDecomposition.PricingSolution(cbdata),
            rc, solvars, solvals, path_data,
        )
    end
    # @show min_rc
    MathOptInterface.submit(
        model.formulation, BlockDecomposition.PricingDualBound(cbdata),
        (stage == 1) ? min_rc : -Inf,
    )
end

# Define the function to separate capacity cuts via RCSP library
function separate_capacity_cuts!(cbdata::CBD, sol::VrpSolution, model::VrpModel) where {CBD}
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
                m.coeff * model.variables_by_id[m.varid] for m in cut.members
            ),
            domain(cut),
        )
        MathOptInterface.submit(
            model.formulation, MathOptInterface.UserCut(cbdata), moi_cut,
        )
    end
    return length(cuts)
end

function separate_rank_one_cuts!(cbdata::CBD, sol::VrpSolution, model::VrpModel) where {CBD}
    println("Separating Rank-1 Cuts with up to $(model.cutsep_phase) rows")
    cuts = run_rank_one_cut_separator(model, sol)

    # add the separated cuts
    for cut in cuts
        moi_cut = ScalarConstraint(
            JuMP.AffExpr(0.0), MathOptInterface.LessThan(cut.rhs),
        )
        MathOptInterface.submit(
            model.formulation, MathOptInterface.UserCut(cbdata), moi_cut, cut,
        )
    end
    return length(cuts)
end

function separate_all_cuts!(cbdata::CBD, model::VrpModel) where {CBD}
    # Check if cuts should be separated
    storage = Coluna.MathProg.getstorage(cbdata.form)
    unit = storage.units[VrpNodeInfoUnit].storage_unit
    if unit.should_stop_cutsep
        println("----- Cut generation is stopped due to tailing off -----")
        unit.should_stop_cutsep = false
        return nothing
    end
    if should_solve_by_mip(unit, model)
        return nothing
    end

    # get the solution of the current relaxation
    sol = VrpSolution(Tuple{Float64, PathVarData}[])
    for (varid, varval) in cbdata.orig_sol
        var = Coluna.MathProg.getvar(cbdata.form, varid)
        if typeof(var.custom_data) == PathVarData
            push!(sol.paths, (varval, var.custom_data))
        end
    end

    # Separate the cuts
    max_rows = get_rcsp_rank1cut_param_value(Int, model.parameters[1], :maxNumRows)
    max_cuts = get_rcsp_rank1cut_param_value(Int, model.parameters[1], :maxNumPerRound)
    if isempty(model.rcc_separators)
        model.cutsep_phase = min(3, max_rows)
    end
    nb_cuts = separate_capacity_cuts!(cbdata, sol, model)
    if (nb_cuts == 0) && (model.cutsep_phase == 0)
        model.cutsep_phase = min(3, max_rows)
    end
    if model.cutsep_phase >= 1 && max_rows > 0 && max_cuts > 0
        nb_new_cuts = separate_rank_one_cuts!(cbdata, sol, model)
        nb_cuts += nb_new_cuts
    end
    if nb_cuts > 0
        unit.separated_cuts = true
    end
    # @show nb_cuts
    return nothing

    # TODO: call the rank one cut separation
    # - ccall `columnGenerationTerminated(false, 0, 0, -1, 0.0, 0.0, false, bool&)`
    #   to check whether to stop rank one cut separation by pricing time
end

function VrpOptimizer(model::VrpModel, config_fname::String, _::AbstractString)
    model.cfg_fname = config_fname
    empty!(model.parameters)
    push!(model.parameters, VrpParameters(config_fname))
    print_params(model.parameters[1].coluna_vrp_params)
    build_capacity_cut_separators!(model)
    max_cuts = get_rcsp_rank1cut_param_value(Int, model.parameters[1], :maxNumPerRound)
    max_rows = get_rcsp_rank1cut_param_value(Int, model.parameters[1], :maxNumRows)
    if max_cuts > 0 && max_rows > 0
        build_rank_one_cut_separator!(model)
    end
    build_solvers!(model)

    # apply the decomposition and store the axis
    nb_graphs = length(model.rcsp_instances)
    @axis(VrpGraphs, 1:nb_graphs)
    @dantzig_wolfe_decomposition(model.formulation, decomp, VrpGraphs)
    model.bd_graphs[1] = decomp

    # Register the custom variable and constraint data
    customvars!(model.formulation, PathVarData)
    customconstrs!(model.formulation, RankOneCutData)

    # Set the mapped variables as representatives
    subproblems = getsubproblems(decomp)
    for var in model.variables_by_id
        subproblemrepresentative(var, subproblems)
    end

    # set the cut callback
    MathOptInterface.set(
        model.formulation, MathOptInterface.UserCutCallback(),
        (cbdata -> separate_all_cuts!(cbdata, model)),
    )

    # preallocate a vector to store the reduced costs
    rcosts = zeros(Float64, get_maxvarid(model))

    # set the solution multiplicities and the pricing callback for each graph
    for g in 1:nb_graphs
        (L, U) = model.rcsp_instances[g].graph.bounds
        specify!(
            subproblems[g], lower_multiplicity = L, upper_multiplicity = U,
            solver = [
                cbdata -> solve_RCSP_pricing(cbdata, stage, model, rcosts) for stage in 1:3
            ],
        )
    end

    # set the callback for reduced-cost fixing and enumeration, and solver by MIP
    model.redcostfix_enum_algo.func =
        (masterform, cbdata, optstate) -> run_redcostfixing_and_enumeration!(
            masterform, cbdata, model, rcosts, optstate,
        )
    model.solve_by_mip_algo.func =
        (masterform, cbdata_vec, optstate) -> run_solve_by_mip!(
            masterform, model, cbdata_vec, optstate,
        )

    # create a dictionary of return values for compatibility with old applications
    sd = Dict(
        MathOptInterface.OPTIMAL => :Optimal,
        MathOptInterface.INFEASIBLE => :Infeasible,
        MathOptInterface.ITERATION_LIMIT => :UserLimit,
        MathOptInterface.TIME_LIMIT => :UserLimit,
        MathOptInterface.OPTIMIZE_NOT_CALLED => :NotSolved,
    )

    # create the optimizer and return it
    return VrpOptimizer(model, sd)
end

function set_cutoff!(opt::VrpOptimizer, cutoffvalue::Float64)
    objectiveprimalbound!(opt.model.formulation, cutoffvalue)
    opt.model.cutoffvalue = cutoffvalue
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
    return JuMP.value(var)
end
