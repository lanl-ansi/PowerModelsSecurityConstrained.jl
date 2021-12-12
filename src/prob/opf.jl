

"""
An OPF formulation conforming to the ARPA-e GOC Challenge 1 specification.
Power balance and branch flow constraints are strictly enforced.
The primary departure from the PowerModels standard formulation is dispatchable
bus shunts and a slight change in the transformer model.
"""
function run_c1_opf_shunt(file, model_constructor, solver; kwargs...)
    return _PM.run_model(file, model_constructor, solver, build_c1_opf_shunt; ref_extensions=[ref_c1!], kwargs...)
end

function build_c1_opf_shunt(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)

    variable_c1_shunt_admittance_imaginary(pm)

    _PM.objective_min_fuel_cost(pm)

    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_c1_power_balance_shunt_dispatch(pm, i)
    end

    for (i,branch) in ref(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end
end


"""
An OPF formulation conforming to the ARPA-e GOC Challenge 1 specification.
Power balance are strictly enforced and the branch flow violations are
penalized based on a conservative linear approximation of the formulation's
flow violation penalty specification.
"""
function run_c1_opf_cheap(file, model_constructor, solver; kwargs...)
    return _PM.run_model(file, model_constructor, solver, build_c1_opf_cheap; ref_extensions=[ref_c1!], kwargs...)
end


function build_c1_opf_cheap(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm, bounded=false)

    variable_c1_branch_power_slack(pm)
    variable_c1_shunt_admittance_imaginary(pm)

    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_c1_power_balance_shunt_dispatch(pm, i)
    end

    for (i,branch) in ref(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        constraint_c1_thermal_limit_from_soft(pm, i)
        constraint_c1_thermal_limit_to_soft(pm, i)
    end

    ##### Setup Objective #####
    _PM.objective_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = var(pm, :pg_cost)
    sm_slack = var(pm, :sm_slack)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) ) +
        sum( 5e5*sm_slack[i] for (i,branch) in ref(pm, :branch_sm_active) )
    )
end


"""
A variant of `run_opf_cheap` model, specialized for solving very large
AC Power Flow models in rectangular coordinates for faster derivative
computations.  Support sparse collections of flow constrains for
increased performance.
"""
function run_c1_opf_cheap_lazy_acr(file, solver; kwargs...)
    return _PM.run_model(file, _PM.ACRPowerModel, solver, build_c1_opf_cheap_lazy_acr; ref_extensions=[ref_c1!], kwargs...)
end

""
function build_c1_opf_cheap_lazy_acr(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm, bounded=false)
    _PM.variable_gen_power(pm)

    variable_c1_branch_power_slack(pm)
    variable_c1_shunt_admittance_imaginary(pm)

    variable_c1_bus_voltage_magnitude_delta(pm)
    variable_c1_gen_power_real_delta(pm)

    vvm = var(pm)[:vvm] = @variable(pm.model,
        [i in ids(pm, :bus)], base_name="vvm",
        lower_bound = ref(pm, :bus, i, "vmin")^2,
        upper_bound = ref(pm, :bus, i, "vmax")^2,
        start = 1.0
    )


    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :branch)
        expression_c1_branch_power_ohms_yt_from(pm, i)
        _PM.expression_branch_power_ohms_yt_to(pm, i)
    end


    vr = var(pm, :vr)
    vi = var(pm, :vi)
    for (i,bus) in ref(pm, :bus)
        vm_midpoint = (bus["vmax"] + bus["vmin"])/2.0
        vm_target = min(vm_midpoint + 0.04, bus["vmax"])
        #vm_target = bus["vm"]

        @constraint(pm.model, vr[i]^2 + vi[i]^2 == vvm[i])
        @constraint(pm.model, vvm[i] == vm_target^2 + var(pm, :vvm_delta, i))
    end

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end
    #Memento.info(_LOGGER, "misc constraints time: $(time() - start_time)")


    start_time = time()
    for (i,branch) in ref(pm, :branch)
        if haskey(branch, "rate_a")
            f_bus_id = branch["f_bus"]
            t_bus_id = branch["t_bus"]
            f_idx = (i, f_bus_id, t_bus_id)
            t_idx = (i, t_bus_id, f_bus_id)

            p_fr = @variable(pm.model, base_name="p_fr", start = 0.0)
            q_fr = @variable(pm.model, base_name="q_fr", start = 0.0)
            p_to = @variable(pm.model, base_name="p_to", start = 0.0)
            q_to = @variable(pm.model, base_name="q_to", start = 0.0)

            @constraint(pm.model, var(pm, :p, f_idx) == p_fr)
            @constraint(pm.model, var(pm, :q, f_idx) == q_fr)
            @constraint(pm.model, var(pm, :p, t_idx) == p_to)
            @constraint(pm.model, var(pm, :q, t_idx) == q_to)

            rating = branch["rate_a"]
            sm_slack = var(pm, :sm_slack, i)
            JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (rating + sm_slack)^2)
            JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (rating + sm_slack)^2)
        end
    end
    #Memento.info(_LOGGER, "flow expr time: $(time() - start_time)")


    start_time = time()
    for (i,gen) in ref(pm, :gen)
        constraint_c1_gen_power_real_deviation(pm, i)
    end
    #Memento.info(_LOGGER, "gen expr time: $(time() - start_time)")


    p = var(pm, :p)
    q = var(pm, :q)
    pg = var(pm, :pg)
    qg = var(pm, :qg)
    bs = var(pm, :bs)
    for (i,bus) in ref(pm, :bus)
        #_PM.constraint_power_balance(pm, i)

        bus_arcs = ref(pm, :bus_arcs, i)
        bus_gens = ref(pm, :bus_gens, i)
        bus_loads = ref(pm, :bus_loads, i)
        bus_shunts_const = ref(pm, :bus_shunts_const, i)
        bus_shunts_var = ref(pm, :bus_shunts_var, i)

        bus_pd = Dict(k => ref(pm, :load, k, "pd") for k in bus_loads)
        bus_qd = Dict(k => ref(pm, :load, k, "qd") for k in bus_loads)

        bus_gs_const = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts_const)
        bus_bs_const = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts_const)

        @constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vvm[i])
        @constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vvm[i] + sum(bs[s]*vvm[i] for s in bus_shunts_var))
    end
    #Memento.info(_LOGGER, "power balance constraint time: $(time() - start_time)")

    vvm_delta = var(pm, :vvm_delta)
    sm_slack = var(pm, :sm_slack)
    pg_delta = var(pm, :pg_delta)

    @objective(pm.model, Min,
        sum( 1e7*vvm_delta[i]^2 for (i,bus) in ref(pm, :bus)) +
        sum( 5e5*sm_slack[i] for (i,branch) in ref(pm, :branch_sm_active)) +
        sum( 1e5*pg_delta[i]^2 for (i,gen) in ref(pm, :gen))
    )
end



""
function run_c1_opf_cheap_target_acp(file, solver; kwargs...)
    return _PM.run_model(file, _PM.ACPPowerModel, solver, build_c1_opf_cheap_target_acp; ref_extensions=[ref_c1!], kwargs...)
end

""
function build_c1_opf_cheap_target_acp(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm, bounded=false)

    variable_c1_branch_power_slack(pm)
    variable_c1_shunt_admittance_imaginary(pm)

    variable_c1_bus_voltage_magnitude_delta(pm)
    variable_c1_gen_power_real_delta(pm)

    _PM.constraint_model_voltage(pm)

    vm = var(pm, :vm)
    for (i,bus) in ref(pm, :bus)
        vm_midpoint = (bus["vmax"] + bus["vmin"])/2.0
        vm_target = min(vm_midpoint + 0.04, bus["vmax"])

        @constraint(pm.model, vm[i] == vm_target + var(pm, :vvm_delta, i))
    end

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    start_time = time()
    for (i,gen) in ref(pm, :gen)
        constraint_c1_gen_power_real_deviation(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_c1_power_balance_shunt_dispatch(pm, i)
    end

    for (i,branch) in ref(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end


    vvm_delta = var(pm, :vvm_delta)
    sm_slack = var(pm, :sm_slack)
    pg_delta = var(pm, :pg_delta)

    @objective(pm.model, Min,
        sum( 1e8*vvm_delta[i]^2 for (i,bus) in ref(pm, :bus)) +
        sum( 5e5*sm_slack[i] for (i,branch) in ref(pm, :branch)) +
        sum( 1e5*pg_delta[i]^2 for (i,gen) in ref(pm, :gen))
    )
end



"""
An OPF formulation conforming to the ARPA-e GOC Challenge 2 specification.

The primary departure from the PowerModels standard formulations are: (1) the
addition of price sensitive loads with a corresponding value maximizing
objective; (2) ramping constraints from a previous operating points; (3) the
addition of penalized slack values that are used for power balance and branch
flow limits, i.e. soft constraints.
"""
function run_c2_opf_soft(file, model_type::Type, optimizer; kwargs...)
    return _PM.run_model(file, model_type, optimizer, build_c2_opf_soft; kwargs...)
end

""
function build_c2_opf_soft(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    variable_bus_delta_abs(pm)
    _PM.variable_gen_power(pm)
    variable_c2_load_power_factor_range(pm)
    _PM.variable_branch_power(pm)
    variable_c2_branch_limit_slack(pm)


    _PM.constraint_model_voltage(pm)
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

    for (i, branch) in ref(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        va_fr = get(ref(pm, :bus, branch["f_bus"]), "va_start", 0.0)
        va_to = get(ref(pm, :bus, branch["t_bus"]), "va_start", 0.0)
        va_detla = va_fr - va_to
        if va_detla <= branch["angmax"] && va_detla >= branch["angmin"]
            _PM.constraint_voltage_angle_difference(pm, i)
        else
            warn(_LOGGER, "skipping constraint_voltage_angle_difference on branch $(i) due to va delta of $(va_detla) for given bounds $(branch["angmin"]) - $(branch["angmax"])")
        end

        constraint_c2_flow_limit_from_soft(pm, i)
        constraint_c2_flow_limit_to_soft(pm, i)
    end


    _PM.objective_variable_pg_cost(pm)
    objective_c2_variable_pd_value(pm)

    p_vio_cost = ref(pm, :p_delta_cost_approx)
    q_vio_cost = ref(pm, :q_delta_cost_approx)
    sm_vio_cost = ref(pm, :sm_cost_approx)

    delta = ref(pm, :delta)

    @objective(pm.model, Max,
        delta * sum(
            sum( var(pm, n, :pd_value, i) for (i,load) in nw_ref[:load])
            - sum( var(pm, n, :pg_cost, i) + gen["oncost"] for (i,gen) in nw_ref[:gen])
            - sum( sm_vio_cost*var(pm, n, :sm_slack, i) for (i,branch) in nw_ref[:branch])
            - sum( p_vio_cost*var(pm, n, :p_delta_abs, i) for (i,bus) in nw_ref[:bus])
            - sum( q_vio_cost*var(pm, n, :q_delta_abs, i) for (i,bus) in nw_ref[:bus])
        for (n, nw_ref) in _PM.nws(pm))
    )

end


"""
An OPF formulation conforming to the ARPA-e GOC Challenge 2 specification.

This formulation is design for optimizing the second-stage contingency problems
It is identical to `run_c2_opf_soft` except for the using different data
parameters which may change between the basecase and contingencies, most
notably the ramping constraints and time passage constant `delta`.
"""
function run_c2_opf_soft_ctg(file, model_type::Type, optimizer; kwargs...)
    return _PM.run_model(file, model_type, optimizer, build_c2_opf_soft_ctg; kwargs...)
end

""
function build_c2_opf_soft_ctg(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    variable_bus_delta_abs(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    variable_c2_branch_limit_slack(pm)
    variable_c2_load_power_factor_range(pm)


    _PM.constraint_model_voltage(pm)
    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    delta_r = ref(pm, :deltarctg)

    for (i, load) in ref(pm, :load)
        z_demand = var(pm, :z_demand, i)

        pd_min = load["pd_nominal"]*JuMP.lower_bound(z_demand)
        pd_max = load["pd_nominal"]*JuMP.upper_bound(z_demand)

        ramping_lb = load["pd_prev"] - load["prdmaxctg"]*delta_r
        ramping_ub = load["pd_prev"] + load["prumaxctg"]*delta_r

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

        ramping_lb = pg_prev - delta_r*gen["prdmaxctg"]
        ramping_ub = pg_prev + delta_r*gen["prumaxctg"]
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

    for (i,branch) in ref(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        #constraint_voltage_angle_difference(pm, i)
        va_fr = get(ref(pm, :bus, branch["f_bus"]), "va_start", 0.0)
        va_to = get(ref(pm, :bus, branch["t_bus"]), "va_start", 0.0)
        va_detla = va_fr - va_to
        if va_detla <= branch["angmax"] && va_detla >= branch["angmin"]
            _PM.constraint_voltage_angle_difference(pm, i)
        else
            warn(_LOGGER, "skipping constraint_voltage_angle_difference on branch $(i) due to va delta of $(va_detla) for given bounds $(branch["angmin"]) - $(branch["angmax"])")
        end

        constraint_c2_flow_limit_from_soft(pm, i)
        constraint_c2_flow_limit_to_soft(pm, i)
    end


    _PM.objective_variable_pg_cost(pm)
    objective_c2_variable_pd_value(pm)

    p_vio_cost = ref(pm, :p_delta_cost_approx)
    q_vio_cost = ref(pm, :q_delta_cost_approx)
    sm_vio_cost = ref(pm, :sm_cost_approx)

    delta = ref(pm, :deltactg)

    @objective(pm.model, Max,
        delta * sum(
            sum( var(pm, n, :pd_value, i) for (i,load) in nw_ref[:load])
            - sum( var(pm, n, :pg_cost, i) + gen["oncost"] for (i,gen) in nw_ref[:gen])
            - sum( sm_vio_cost*var(pm, n, :sm_slack, i) for (i,branch) in nw_ref[:branch])
            - sum( p_vio_cost*var(pm, n, :p_delta_abs, i) for (i,bus) in nw_ref[:bus])
            - sum( q_vio_cost*var(pm, n, :q_delta_abs, i) for (i,bus) in nw_ref[:bus])
        for (n, nw_ref) in _PM.nws(pm))
    )
end



"""
An OPF formulation with generator unit commitment conforming to the ARPA-e GOC
Challenge 2 specification.

The primary departure from the PowerModels standard formulations are: (1) the
addition of price sensitive loads with a corresponding value maximizing
objective; (2) ramping constraints from a previous operating points; (3) the
addition of penalized slack values for branch flow limits, i.e. soft
constraints.

This model differs from the `c2_opf_soft` model by strictly enforcing power
balance at the network buses.
"""
function run_c2_opf_uc(file, model_type::Type, optimizer; kwargs...)
    return _PM.run_model(file, model_type, optimizer, build_c2_opf_uc; kwargs...)
end

""
function build_c2_opf_uc(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)

    _PM.variable_gen_indicator(pm)
    _PM.variable_gen_power_on_off(pm)
    nw=nw_id_default
    gen_su = var(pm, nw)[:gen_su] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_gen_su",
        binary = true,
        start = _PM.comp_start_value(ref(pm, nw, :gen, i), "gen_su_start", 0.0)
    )
    _PM.sol_component_value(pm, nw, :gen, :gen_su, ids(pm, nw, :gen), gen_su)

    gen_sd = var(pm, nw)[:gen_sd] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_gen_sd",
        binary = true,
        start = _PM.comp_start_value(ref(pm, nw, :gen, i), "gen_sd_start", 0.0)
    )
    _PM.sol_component_value(pm, nw, :gen, :gen_sd, ids(pm, nw, :gen), gen_sd)


    variable_c2_load_power_factor_range(pm)
    _PM.variable_branch_power(pm)
    variable_c2_branch_limit_slack(pm)


    _PM.constraint_model_voltage(pm)
    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    delta_r = ref(pm, :deltar)

    for (i, load) in ref(pm, :load)
        z_demand = var(pm, :z_demand, i)
        
        @constraint(pm.model, load["pd_nominal"]*z_demand <= load["pd_prev"] + load["prumax"]*delta_r)
        @constraint(pm.model, load["pd_nominal"]*z_demand >= load["pd_prev"] - load["prdmax"]*delta_r)
    end

    for (i, gen) in ref(pm, :gen)
        _PM.constraint_gen_power_on_off(pm, i)

        pg = var(pm, :pg, i)
        pg_prev = var(pm, :z_gen, i)*gen["pg_prev"]
        if gen["status_prev"] == 0 # starting up
           pg_prev = var(pm, :z_gen, i)*gen["pmin"]
        end

        @constraint(pm.model, pg <= pg_prev + delta_r*gen["prumax"])
        @constraint(pm.model, pg >= pg_prev - delta_r*gen["prdmax"])

        @constraint(pm.model, var(pm, :z_gen, i) - gen["status_prev"] == var(pm, :gen_su, i) - var(pm, :gen_sd, i))
        @constraint(pm.model, var(pm, :gen_su, i) + var(pm, :gen_sd, i) <= 1)

        if gen["status_prev"] == 0 && gen["suqual"] == 0
            @constraint(pm.model, var(pm, :z_gen, i) == 0)
        end

        if gen["status_prev"] == 1 && gen["sdqual"] == 0
            @constraint(pm.model, var(pm, :z_gen, i) == 1)
        end

        # if haskey(gen, "gen_status_fixed")
        #     @constraint(pm.model, var(pm, :z_gen, i) == gen["gen_status_fixed"])
        # end
    end

    for i in ids(pm, :bus)
        constraint_c2_power_balance(pm, i)
    end

    for (i, branch) in ref(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        #constraint_voltage_angle_difference(pm, i)
        va_fr = get(ref(pm, :bus, branch["f_bus"]), "va_start", 0.0)
        va_to = get(ref(pm, :bus, branch["t_bus"]), "va_start", 0.0)
        va_detla = va_fr - va_to
        if va_detla <= branch["angmax"] && va_detla >= branch["angmin"]
            _PM.constraint_voltage_angle_difference(pm, i)
        else
            warn(_LOGGER, "skipping constraint_voltage_angle_difference on branch $(i) due to va delta of $(va_detla) for given bounds $(branch["angmin"]) - $(branch["angmax"])")
        end

        constraint_c2_flow_limit_from_soft(pm, i)
        constraint_c2_flow_limit_to_soft(pm, i)
    end


    _PM.objective_variable_pg_cost(pm)
    objective_c2_variable_pd_value(pm)

    #p_vio_cost = ref(pm, :p_delta_cost_approx)
    #q_vio_cost = ref(pm, :q_delta_cost_approx)
    sm_vio_cost = ref(pm, :sm_cost_approx)

    delta = ref(pm, :delta)

    gen_on_cost = Dict(i => 1.0*gen["oncost"] for (i, gen) in ref(pm, :gen))

    @objective(pm.model, Max,
        delta * sum(
            sum( var(pm, n, :pd_value, i) for (i,load) in nw_ref[:load])
            - sum( 
                var(pm, n, :pg_cost, i) + gen_on_cost[i]*var(pm, :z_gen, i) +
                (gen["sucost"]/delta)*var(pm, n, :gen_su, i) + 
                (gen["sdcost"]/delta)*var(pm, n, :gen_sd, i) 
            for (i,gen) in nw_ref[:gen])
            - sum( sm_vio_cost*var(pm, n, :sm_slack, i) for (i,branch) in nw_ref[:branch])
            #- sum( p_vio_cost*var(pm, n, :p_delta_abs, i) for (i,bus) in nw_ref[:bus])
            #- sum( q_vio_cost*var(pm, n, :q_delta_abs, i) for (i,bus) in nw_ref[:bus])
        for (n, nw_ref) in _PM.nws(pm))
    )
end





