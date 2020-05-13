
"Solves DC-SCOPF by adding branch-based PTDF cuts to the base-case iterativly"
function run_scopf_ptdf_cuts(file::String, model_type::Type, optimizer; kwargs...)
    data = _PM.parse_file(file)
    return run_scopf_ptdf_cuts!(data, optimizer; kwargs...)
end

function run_scopf_ptdf_cuts!(network::Dict{String,<:Any}, model_type::Type, optimizer; max_iter::Int=100, time_limit::Float64=Inf)
    time_start = time()

    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = _PM.run_opf(network, model_type, optimizer)
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        error(_LOGGER, "base-case OPF solve failed in run_scopf_ptdf_cuts, status $(result["termination_status"])")
    end
    _PM.update_data!(network, result["solution"])

    result["iterations"] = 0

    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = check_contingencies_branch_power(network, total_cut_limit=iteration, gen_flow_cuts=[], branch_flow_cuts=[])

        cuts_found = length(cuts.gen_cuts)+length(cuts.branch_cuts)
        if cuts_found <= 0
            info(_LOGGER, "no violated cuts found scopf fixed-point reached")
            break
        else
            info(_LOGGER, "found $(length(cuts.gen_cuts) + length(cuts.branch_cuts)) branch flow violations")
        end

        append!(network["gen_flow_cuts"], cuts.gen_cuts)
        append!(network["branch_flow_cuts"], cuts.branch_cuts)

        info(_LOGGER, "active cuts: gen $(length(network["gen_flow_cuts"])), branch $(length(network["branch_flow_cuts"]))")

        time_solve_start = time()
        result = run_scopf_cuts(network, model_type, optimizer)
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

