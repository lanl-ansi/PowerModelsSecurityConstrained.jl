function run_c1_scopf_contigency_cuts(ini_file::String, model_type::Type, optimizer; scenario_id::String="", kwargs...)
    goc_data = parse_c1_case(ini_file, scenario_id=scenario_id)
    network = build_c1_pm_model(goc_data)
    return run_c1_scopf_contigency_cuts(network, model_type, optimizer; kwargs...)
end

"""
Solves a SCOPF problem by iteratively checking for violated contingencies and
resolving until a fixed-point is reached
"""
function run_c1_scopf_contigency_cuts(network::Dict{String,<:Any}, model_type::Type, optimizer; max_iter::Int=100, time_limit::Float64=Inf)
    if _IM.ismultinetwork(network)
        error(_LOGGER, "run_c1_scopf_contigency_cuts can only be used on single networks")
    end

    time_start = time()

    network_base = deepcopy(network)
    network_active = deepcopy(network)

    gen_contingencies = network_base["gen_contingencies"]
    branch_contingencies = network_base["branch_contingencies"]

    network_active["gen_contingencies"] = []
    network_active["branch_contingencies"] = []

    multinetwork = build_c1_scopf_multinetwork(network_active)

    result = run_c1_scopf(multinetwork, model_type, optimizer)
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        error(_LOGGER, "base-case SCOPF solve failed in run_c1_scopf_contigency_cuts, status $(result["termination_status"])")
    end
    #_PM.print_summary(result["solution"])
    solution = result["solution"]["nw"]["0"]
    solution["per_unit"] = result["solution"]["per_unit"]

    _PM.update_data!(network_base, solution)
    _PM.update_data!(network_active, solution)

    result["iterations"] = 0

    iteration = 1
    contingencies_found = 1
    while contingencies_found > 0
        time_start_iteration = time()

        contingencies = check_c1_contingency_violations(network_base, contingency_limit=iteration)
        #println(contingencies)

        contingencies_found = 0
        #append!(network_active["gen_contingencies"], contingencies.gen_contingencies)
        for cont in contingencies.gen_contingencies
            if cont in network_active["gen_contingencies"]
                warn(_LOGGER, "generator contingency $(cont.label) is active but not secure")
            else
                push!(network_active["gen_contingencies"], cont)
                contingencies_found += 1
            end
        end

        #append!(network_active["branch_contingencies"], contingencies.branch_contingencies)
        for cont in contingencies.branch_contingencies
            if cont in network_active["branch_contingencies"]
                warn(_LOGGER, "branch contingency $(cont.label) is active but not secure")
            else
                push!(network_active["branch_contingencies"], cont)
                contingencies_found += 1
            end
        end

        if contingencies_found <= 0
            info(_LOGGER, "no new violated contingencies found, scopf fixed-point reached")
            break
        else
            info(_LOGGER, "found $(contingencies_found) new contingencies with violations")
        end


        info(_LOGGER, "active contingencies: gen $(length(network_active["gen_contingencies"])), branch $(length(network_active["branch_contingencies"]))")

        time_solve_start = time()
        multinetwork = build_c1_scopf_multinetwork(network_active)
        result = run_c1_scopf(multinetwork, model_type, optimizer)
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            warn(_LOGGER, "scopf solve failed with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        # for (nw,nw_sol) in result["solution"]["nw"]
        #     if nw != "0"
        #         println(nw, " ", nw_sol["delta"])
        #     end
        # end
        info(_LOGGER, "objective: $(result["objective"])")
        solution = result["solution"]["nw"]["0"]
        solution["per_unit"] = result["solution"]["per_unit"]

        _PM.update_data!(network_base, solution)
        _PM.update_data!(network_active, solution)

        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        iteration += 1
    end

    result["solution"] = solution
    result["iterations"] = iteration

    return result
end


function run_c1_scopf_ptdf_cuts(ini_file::String, model_type::Type, optimizer; scenario_id::String="", kwargs...)
    goc_data = parse_c1_case(ini_file, scenario_id=scenario_id)
    network = build_c1_pm_model(goc_data)
    return run_c1_scopf_ptdf_cuts!(network, model_type, optimizer; kwargs...)
end

"""
Solves a SCOPF problem by iteratively checking for violated branch flow
constraints in contingencies and resolving until a fixed-point is reached.

The base-case model is formulation agnostic.  The flow cuts are based on PTDF
and utilize the DC Power Flow assumption.
"""
function run_c1_scopf_ptdf_cuts!(network::Dict{String,<:Any}, model_type::Type, optimizer; max_iter::Int=100, time_limit::Float64=Inf)
    if _IM.ismultinetwork(network)
        error(_LOGGER, "run_c1_scopf_ptdf_cuts can only be used on single networks")
    end

    time_start = time()

    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = _PM.solve_opf(network, model_type, optimizer)
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        error(_LOGGER, "base-case OPF solve failed in run_c1_scopf_ptdf_cuts, status $(result["termination_status"])")
    end
    info(_LOGGER, "objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"])

    result["iterations"] = 0

    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = check_c1_contingencies_branch_power(network, total_cut_limit=iteration, gen_flow_cuts=[], branch_flow_cuts=[])

        cuts_found = length(cuts.gen_cuts)+length(cuts.branch_cuts)
        if cuts_found <= 0
            info(_LOGGER, "no violated cuts found scopf fixed-point reached")
            break
        else
            info(_LOGGER, "found $(cuts_found) branch flow violations")
        end

        append!(network["gen_flow_cuts"], cuts.gen_cuts)
        append!(network["branch_flow_cuts"], cuts.branch_cuts)

        info(_LOGGER, "active cuts: gen $(length(network["gen_flow_cuts"])), branch $(length(network["branch_flow_cuts"]))")

        time_solve_start = time()
        result = run_c1_scopf_cuts(network, model_type, optimizer)
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            warn(_LOGGER, "scopf solve failed with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        info(_LOGGER, "objective: $(result["objective"])")
        _PM.update_data!(network, result["solution"])

        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        iteration += 1
    end

    result["iterations"] = iteration
    return result
end

