
"""
An OTS formulation with topology optimization conforming to the ARPA-e GOC
Challenge 2 specification.

This model is similar in spirit to `c2_opf_soft`, however, branch flow limits
are strictly enforced to improve the reliability of solving the non-convex
optimization problems.
"""
function run_c2_ots_soft_bus(file, model_type::Type, optimizer; kwargs...)
    return _PM.run_model(file, model_type, optimizer, build_c2_ots_soft_bus; ref_extensions=[_PM.ref_add_on_off_va_bounds!], kwargs...)
end

""
function build_c2_ots_soft_bus(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage_on_off(pm)
    variable_bus_delta_abs(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_indicator(pm)
    _PM.variable_branch_power(pm)
    variable_c2_load_power_factor_range(pm)

    nw=nw_id_default
    z_branch_delta_abs = var(pm, nw)[:z_branch_delta_abs] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :branch)], base_name="$(nw)_z_branch_delta_abs",
        binary = true,
        start = _PM.comp_start_value(ref(pm, nw, :branch, i), "z_branch_delta_abs_start", 0.0)
    )
    _PM.sol_component_value(pm, nw, :branch, :z_branch_delta_abs, ids(pm, nw, :branch), z_branch_delta_abs)

    for (i, branch) in ref(pm, :branch)
        if branch["swqual"] == 0
            @constraint(pm.model, var(pm, :z_branch, i) == branch["status_prev"])
            @constraint(pm.model, var(pm, :z_branch_delta_abs, i) == 0)
            #var(pm, :z_branch)[i] = branch["status_prev"]
        else
            @constraint(pm.model, var(pm, :z_branch_delta_abs, i) >=  var(pm, :z_branch, i) - branch["status_prev"])
            @constraint(pm.model, var(pm, :z_branch_delta_abs, i) >= -var(pm, :z_branch, i) + branch["status_prev"])
        end
    end
    var(pm)[:z_branch] = Dict(branch["swqual"] == 0 ? i => branch["status_prev"] : i => var(pm, :z_branch, i) for (i, branch) in ref(pm, :branch))



    _PM.constraint_model_voltage_on_off(pm)
    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    delta_r = ref(pm, :deltar)

    for (i, load) in ref(pm, :load)
        z_demand = var(pm, :z_demand, i)
        
        pd_min = load["pd_nominal"]*JuMP.lower_bound(z_demand)
        pd_max = load["pd_nominal"]*JuMP.upper_bound(z_demand)

        ramping_lb = load["pd_prev"] - load["prdmax"]*delta_r
        ramping_ub = load["pd_prev"] + load["prumax"]*delta_r

        if pd_min < ramping_lb
            #println("update lower bound $(pd_min) -> $(ramping_lb)")
            JuMP.set_lower_bound(z_demand, ramping_lb/load["pd_nominal"])
        end

        if pd_max > ramping_ub
            #println("update upper bound $(pd_max) -> $(ramping_ub)")
            JuMP.set_upper_bound(z_demand, ramping_ub/load["pd_nominal"])
        end
    end

    for (i, gen) in ref(pm, :gen)
        pg = var(pm, :pg, i)

        pg_prev = gen["pg_prev"]
        if gen["status_prev"] == 0 # starting up
           pg_prev = gen["pmin"]
           warn(_LOGGER, "generator $(i) startung up setting previous value to $(pg_prev)")
        end

        ramping_lb = pg_prev - delta_r*gen["prdmax"]
        ramping_ub = pg_prev + delta_r*gen["prumax"]
        if JuMP.lower_bound(pg) < ramping_lb
            #println("update lower bound $(lower_bound(pg)) -> $(ramping_lb)")
            JuMP.set_lower_bound(pg, ramping_lb)
        end

        if JuMP.upper_bound(pg) > ramping_ub
            #println("update upper bound $(upper_bound(pg)) -> $(ramping_ub)")
            JuMP.set_upper_bound(pg, ramping_ub)
        end
    end


    for i in ids(pm, :bus)
        constraint_c2_power_balance_soft_lin(pm, i)
    end

    # can be used to limit the amount of switching considered
    #@constraint(pm.model, sum(var(pm, :z_branch_delta_abs, i) for (i, branch) in ref(pm, :branch)) <= ceil(0.99*length(ref(pm, :branch))))

    for (i, branch) in ref(pm, :branch)
        #constraint_ohms_yt_from_goc(pm, i)
        #constraint_ohms_yt_to(pm, i)
        _PM.constraint_ohms_yt_from_on_off(pm, i)
        _PM.constraint_ohms_yt_to_on_off(pm, i)

        va_fr = get(ref(pm, :bus, branch["f_bus"]), "va_start", 0.0)
        va_to = get(ref(pm, :bus, branch["t_bus"]), "va_start", 0.0)
        va_detla = va_fr - va_to
        if va_detla <= branch["angmax"] && va_detla >= branch["angmin"]
            _PM.constraint_voltage_angle_difference_on_off(pm, i)
        else
            warn(_LOGGER, "skipping constraint_voltage_angle_difference_on_off on branch $(i) due to va delta of $(va_detla) for given bounds $(branch["angmin"]) - $(branch["angmax"])")
        end


        _PM.constraint_thermal_limit_from_on_off(pm, i)
        _PM.constraint_thermal_limit_to_on_off(pm, i)
    end

    _PM.objective_variable_pg_cost(pm)
    objective_c2_variable_pd_value(pm)

    p_vio_cost = ref(pm, :p_delta_cost_approx)
    q_vio_cost = ref(pm, :q_delta_cost_approx)
    #sm_vio_cost = ref(pm, :sm_cost_approx)

    delta = ref(pm, :delta)

    @objective(pm.model, Max,
        delta * sum(
            sum( var(pm, n, :pd_value, i) for (i,load) in nw_ref[:load])
            - sum( var(pm, n, :pg_cost, i) + gen["oncost"] for (i,gen) in nw_ref[:gen])
            - sum( branch["csw"]/delta*var(pm, n, :z_branch_delta_abs, i) for (i,branch) in nw_ref[:branch])
            #- sum( sm_vio_cost*var(pm, n, :sm_slack, i) for (i,branch) in nw_ref[:branch])
            - sum( p_vio_cost*var(pm, n, :p_delta_abs, i) for (i,bus) in nw_ref[:bus])
            - sum( q_vio_cost*var(pm, n, :q_delta_abs, i) for (i,bus) in nw_ref[:bus])
        for (n, nw_ref) in _PM.nws(pm))
    )
end
