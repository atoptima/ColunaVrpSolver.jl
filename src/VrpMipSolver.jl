function solve_vrp_by_mip(
    masterform::Coluna.MathProg.Formulation, model::VrpModel, cbdata_vec::Vector{CB},
    primal_bnd::Float64,
) where {CB}
    mip = Model(
        optimizer_with_attributes(
            CPLEX.Optimizer,
            "MIPGap" => model.parameters[1].coluna_vrp_params.relOptimalityGapTolerance,
            "Cutoff" => primal_bnd,
        ),
    )

    # retrieve the paths and the corresponding mapped variables from all RCSP subproblems
    paths_arcs = Vector{Vector{Int}}[]
    paths_vars = Vector{Dict{VariableRef, Int}}[]
    for rcsp in model.rcsp_instances
        rcsp_paths_arcs = get_enumerated_paths(rcsp)
        push!(paths_arcs, rcsp_paths_arcs)
        rcsp_paths_vars = Dict{VariableRef, Int}[]
        for p in rcsp_paths_arcs
            p_vars = Dict{VariableRef, Int}()
            for a in p
                for var in rcsp.graph.mappings[a+1]
                    p_vars[var] = get(p_vars, var, 0) + 1
                end
            end
            push!(rcsp_paths_vars, p_vars)
        end
        push!(paths_vars, rcsp_paths_vars)
    end

    # declare the paths variables
    nb_rcsps = length(model.rcsp_instances)
    @variable(mip, λ[i in 1:nb_rcsps, j in 1:length(paths_arcs[i])] >= 0, Int)

    # set the objective function
    orig_form_obj = objective_function(model.formulation)
    pathcost(i, j) = sum(
        count * coefficient(orig_form_obj, var) for (var, count) in paths_vars[i][j]
    )
    @objective(mip, Min, sum(
        pathcost(i, j) * λ[i, j] for i in 1:nb_rcsps for j in 1:length(paths_arcs[i])
    ))

    # obtain the master constraints coefficients over the paths
    matrix = Coluna.MathProg.getcoefmatrix(masterform)
    opt = JuMP.unsafe_backend(model.formulation)
    cid_to_coeffs = Dict{Coluna.MathProg.ConstrId, Vector{Tuple{Float64, Int, Int}}}()
    cid_to_sense = Dict{Coluna.MathProg.ConstrId, Symbol}()
    cid_to_rhs = Dict{Coluna.MathProg.ConstrId, Float64}()
    for i in 1:nb_rcsps
        for j in 1:length(paths_arcs[i])
            for (var, count) in paths_vars[i][j]
                vid = Coluna._get_varid_of_origvar_in_form(opt.env, masterform, JuMP.index(var))
                # var_name = Coluna.MathProg.getname(masterform, vid)
                # @show var_name
                for (cid, coeff) in @view matrix[:, vid]
                    if Coluna.MathProg.iscuractive(masterform, cid) &&
                       Coluna.MathProg.getduty(cid) <= Coluna.MathProg.AbstractMasterConstr
                        # ctr_name = Coluna.MathProg.getname(masterform, cid)
                        # @show coeff, ctr_name
                        if !haskey(cid_to_coeffs, cid)
                            coeffs = Tuple{Float64, Int, Int}[]
                            cid_to_coeffs[cid] = coeffs
                            cid_to_rhs[cid] = Coluna.MathProg.getcurrhs(masterform, cid)
                            sense = Coluna.MathProg.getcursense(masterform, cid)
                            if sense == Coluna.MathProg.Equal
                                cid_to_sense[cid] = :(==)
                            elseif sense == Coluna.MathProg.Greater
                                cid_to_sense[cid] = :(>=)
                            else
                                cid_to_sense[cid] = :(<=)
                            end
                        else
                            coeffs = cid_to_coeffs[cid]
                        end
                        push!(coeffs, (count * coeff, i, j))
                    end
                end
            end
        end
    end

    # add the master constraints
    k_to_cid = collect(keys(cid_to_coeffs))
    for k in 1:length(cid_to_coeffs)
        sense = cid_to_sense[k_to_cid[k]]
        if sense == :(==)
            @constraint(mip,
                sum(coeff * λ[i, j] for (coeff, i, j) in cid_to_coeffs[k_to_cid[k]]) ==
                cid_to_rhs[k_to_cid[k]]
            )
        elseif sense == :(>=)
            @constraint(mip,
                sum(coeff * λ[i, j] for (coeff, i, j) in cid_to_coeffs[k_to_cid[k]]) >=
                cid_to_rhs[k_to_cid[k]]
            )
        else
            @constraint(mip,
                sum(coeff * λ[i, j] for (coeff, i, j) in cid_to_coeffs[k_to_cid[k]]) <=
                cid_to_rhs[k_to_cid[k]]
            )
        end
    end

    # add the convexity constraints
    @constraint(mip, conv_lb[i in 1:nb_rcsps],
        sum(λ[i, j] for j in 1:length(paths_arcs[i])) >= model.rcsp_instances[i].graph.bounds[1]
    )
    @constraint(mip, conv_ub[i in 1:nb_rcsps],
        sum(λ[i, j] for j in 1:length(paths_arcs[i])) <= model.rcsp_instances[i].graph.bounds[2]
    )

    # solve the problem and return the solution
    optimize!(mip)
    if JuMP.has_values(mip)
        varids = Coluna.MathProg.VarId[]
        varcoeffs = Float64[]
        cost = 0.0
        for i in 1:nb_rcsps
            for j in 1:length(paths_arcs[i])
                val = value(λ[i, j])
                if val > 0.5
                    # build the subproblem solution data
                    subvarids = [
                        Coluna._get_varid_of_origvar_in_form(
                            opt.env, cbdata_vec[i].form, JuMP.index(var),
                        )
                        for (var, _) in paths_vars[i][j]
                    ]
                    subvarcoeffs = [Float64(count) for (_, count) in paths_vars[i][j]]
                    push!(subvarids, cbdata_vec[i].form.duty_data.setup_var)
                    push!(subvarcoeffs, 1.0)
                    subcost = pathcost(i, j)
                    # subvarnames = [
                    #     Coluna.MathProg.getname(cbdata_vec[i].form, vid) for vid in subvarids
                    # ]
                    # @show subvarnames

                    # add the subproblem solution to the subproblem
                    subsol = Coluna.MathProg.PrimalSolution(
                        cbdata_vec[i].form, subvarids, subvarcoeffs, subcost,
                        Coluna.FEASIBLE_SOL,
                    )
                    col_id = Coluna.MathProg.insert_column!(masterform, subsol, "MC")
                    mc_var = Coluna.MathProg.getvar(masterform, col_id)

                    # accumulate the master solution data
                    push!(varids, Coluna.MathProg.getid(mc_var))
                    push!(varcoeffs, round(val))
                    cost += round(val) * subcost
                end
            end
        end

        # add the master solution to the master problem
        sol = Coluna.MathProg.PrimalSolution(
            masterform, varids, varcoeffs, cost, Coluna.FEASIBLE_SOL,
        )
        return [sol]
    else
        return Coluna.MathProg.PrimalSolution[]
    end
end
