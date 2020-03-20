

"""
An OPF formulation conforming to the ARPA-e GOC Challenge 1 specification.
Power balance and branch flow constraints are strictly enforced.
The primary departure from the PowerModels standard formulation is dispatchable
bus shunts and a slight change in the transformer model.
"""
function run_opf_shunt(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, build_opf_shunt; ref_extensions=[ref_add_goc!], kwargs...)
end

function build_opf_shunt(pm::AbstractPowerModel)
    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm)

    variable_reactive_shunt(pm)

    PowerModels.objective_min_fuel_cost(pm)

    PowerModels.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance_shunt_dispatch(pm, i)
    end

    for (i,branch) in ref(pm, :branch)
        constraint_ohms_yt_from_goc(pm, i)
        PowerModels.constraint_ohms_yt_to(pm, i)

        PowerModels.constraint_voltage_angle_difference(pm, i)

        PowerModels.constraint_thermal_limit_from(pm, i)
        PowerModels.constraint_thermal_limit_to(pm, i)
    end
end


"""
An OPF formulation conforming to the ARPA-e GOC Challenge 1 specification.
Power balance are strictly enforced and the branch flow violations are
penalized based on a conservative linear approximation of the formulation's
flow violation penalty specification.
"""
function run_opf_cheap(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, build_opf_cheap; ref_extensions=[ref_add_goc!], kwargs...)
end


function build_opf_cheap(pm::AbstractPowerModel)
    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm, bounded=false)

    variable_branch_flow_slack(pm)
    variable_reactive_shunt(pm)

    PowerModels.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance_shunt_dispatch(pm, i)
    end

    for (i,branch) in ref(pm, :branch)
        constraint_ohms_yt_from_goc(pm, i)
        PowerModels.constraint_ohms_yt_to(pm, i)

        PowerModels.constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from_soft(pm, i)
        constraint_thermal_limit_to_soft(pm, i)
    end

    ##### Setup Objective #####
    objective_variable_pg_cost(pm)
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
function run_opf_cheap_lazy_acr(file, solver; kwargs...)
    return run_model(file, ACRPowerModel, solver, build_opf_cheap_lazy_acr; ref_extensions=[ref_add_goc!], kwargs...)
end

""
function build_opf_cheap_lazy_acr(pm::AbstractPowerModel)
    PowerModels.variable_voltage(pm, bounded=false)
    PowerModels.variable_generation(pm)

    variable_branch_flow_slack(pm)
    variable_reactive_shunt(pm)

    variable_vvm_delta(pm)
    variable_pg_delta(pm)

    vvm = var(pm)[:vvm] = @variable(pm.model,
        [i in ids(pm, :bus)], base_name="vvm",
        lower_bound = ref(pm, :bus, i, "vmin")^2,
        upper_bound = ref(pm, :bus, i, "vmax")^2,
        start = 1.0
    )


    PowerModels.constraint_model_voltage(pm)

    for i in ids(pm, :branch)
        expression_branch_power_yt_from_goc(pm, i)
        expression_branch_power_yt_to(pm, i)
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
        PowerModels.constraint_theta_ref(pm, i)
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
        constraint_gen_active_deviation(pm, i)
    end
    #Memento.info(_LOGGER, "gen expr time: $(time() - start_time)")


    p = var(pm, :p)
    q = var(pm, :q)
    pg = var(pm, :pg)
    qg = var(pm, :qg)
    bs = var(pm, :bs)
    for (i,bus) in ref(pm, :bus)
        #PowerModels.constraint_power_balance(pm, i)

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
