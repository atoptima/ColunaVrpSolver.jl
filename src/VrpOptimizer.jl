mutable struct VrpOptimizer
    model::VrpModel
    status_dict::Dict{MathOptInterface.TerminationStatusCode,Symbol}
end

function VrpOptimizer(model::VrpModel, config_fname::String, _::AbstractString="")
    model.cfg_fname = config_fname

    # apply the decomposition and store the axis
    nb_graphs = length(model.rcsp_instances)
    @axis(VrpGraphs, 1:nb_graphs)
    @dantzig_wolfe_decomposition(model.formulation, decomp, VrpGraphs)
    model.bd_graphs[1] = decomp

    # Set the mapped variables as representatives
    subproblems = getsubproblems(decomp)
    for var in model.variables_by_id
        subproblemrepresentative(var, subproblems)
    end

    # preallocate a vector to store the reduced costs
    rcosts = zeros(Float64, get_maxvarid(model))

    # set the solution multiplicities and the pricing callback for each graph
    for g in 1:nb_graphs
        (L, U) = model.rcsp_instances[g].graph.bounds
        specify!(
            subproblems[g], lower_multiplicity=L, upper_multiplicity=U,
        )
    end

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
