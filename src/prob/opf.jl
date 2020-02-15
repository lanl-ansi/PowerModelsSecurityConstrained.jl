function solution_goc!(pm::GenericPowerModel, sol::Dict{String,Any})
    PowerModels.add_setpoint_bus_voltage!(sol, pm)
    PowerModels.add_setpoint_generator_power!(sol, pm)
    PowerModels.add_setpoint_branch_flow!(sol, pm)

    PowerModels.add_setpoint!(sol, pm, "branch", "sm_slack", :sm_slack, status_name="br_status")
    PowerModels.add_setpoint!(sol, pm, "bus", "p_delta_abs", :p_delta_abs, status_name="bus_type", inactive_status_value = 4)
    PowerModels.add_setpoint!(sol, pm, "bus", "q_delta_abs", :q_delta_abs, status_name="bus_type", inactive_status_value = 4)

    PowerModels.add_setpoint!(sol, pm, "bus", "vm_offset", :vm_offset, status_name="bus_type", inactive_status_value = 4)

    #PowerModels.add_setpoint!(sol, pm, "gen", "pg_cost", :pg_cost, status_name="gen_status", inactive_status_value = 0, conductorless=true)

    #PowerModels.add_setpoint!(sol, pm, "shunt", "gs", :qsh)
    # requires check dispatchable flag
    #PowerModels.add_setpoint!(sol, pm, "shunt", "bs", :bsh)

    add_setpoint_dispatchable(sol, pm, "shunt", "bs", :bsh, default_value = (item) -> item["bs"], dispatchable_check=true)
end


function run_opf_shunt(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_opf_shunt; ref_extensions=[ref_add_goc!], solution_builder = solution_goc!, kwargs...)
end

function post_opf_shunt(pm::GenericPowerModel)
    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm)
    PowerModels.variable_dcline_flow(pm)

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

    for i in ids(pm, :dcline)
        PowerModels.constraint_dcline(pm, i)
    end
end


function run_opf_cheap(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_opf_cheap; ref_extensions=[ref_add_goc!], solution_builder = solution_goc!, kwargs...)
end


function post_opf_cheap(pm::GenericPowerModel)
    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm, bounded=false)
    PowerModels.variable_dcline_flow(pm)

    variable_reactive_shunt(pm)

    sm_slack = var(pm)[:sm_slack] = @variable(pm.model,
        [l in ids(pm, :branch)], base_name="sm_slack",
        lower_bound = 0.0,
        start = 0.0
    )

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

        #PowerModels.constraint_thermal_limit_from(pm, i)
        #PowerModels.constraint_thermal_limit_to(pm, i)

        f_bus_id = branch["f_bus"]
        t_bus_id = branch["t_bus"]
        f_idx = (i, f_bus_id, t_bus_id)
        t_idx = (i, t_bus_id, f_bus_id)

        rate_a = branch["rate_a"]
        @constraint(pm.model, var(pm, :p, f_idx)^2 + var(pm, :q, f_idx)^2 <= (rate_a + sm_slack[i])^2)
        @constraint(pm.model, var(pm, :p, t_idx)^2 + var(pm, :q, t_idx)^2 <= (rate_a + sm_slack[i])^2)
    end

    ##### Setup Objective #####
    objective_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = var(pm, pm.cnw, :pg_cost)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) ) +
        sum( 5e5*sm_slack[i] for (i,branch) in ref(pm, :branch) )
    )
end


function run_opf_cheap_dc(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_opf_cheap_dc; ref_extensions=[ref_add_goc!], solution_builder = solution_goc!, kwargs...)
end


function post_opf_cheap_dc(pm::GenericPowerModel)
    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm, bounded=false)
    PowerModels.variable_dcline_flow(pm)

    sm_slack = var(pm)[:sm_slack] = @variable(pm.model,
        [l in ids(pm, :branch)], base_name="sm_slack",
        lower_bound = 0.0,
        start = 0.0
    )

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

        #PowerModels.constraint_thermal_limit_from(pm, i)
        #PowerModels.constraint_thermal_limit_to(pm, i)

        f_bus_id = branch["f_bus"]
        t_bus_id = branch["t_bus"]
        f_idx = (i, f_bus_id, t_bus_id)
        t_idx = (i, t_bus_id, f_bus_id)

        rate_a = branch["rate_a"]
        @constraint(pm.model, var(pm, :p, f_idx) <= rate_a + sm_slack[i])
        @constraint(pm.model, var(pm, :p, f_idx) >= -(rate_a + sm_slack[i]))
    end

    ##### Setup Objective #####
    objective_variable_pg_cost(pm)

    # explicit network id needed because of conductor-less
    pg_cost = var(pm, pm.cnw, :pg_cost)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) ) +
        sum( 5e5*sm_slack[i] for (i,branch) in ref(pm, :branch) )
    )
end



""
function run_opf_pg_pf_rect_5(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_opf_pg_pf_rect_5; ref_extensions=[ref_add_goc!], solution_builder=solution_goc!, kwargs...)
end

""
function post_opf_pg_pf_rect_5(pm::GenericPowerModel)
    PowerModels.variable_voltage(pm, bounded=false)
    PowerModels.variable_generation(pm)

    variable_reactive_shunt(pm)

    vvm = var(pm)[:vvm] = @variable(pm.model,
        [i in ids(pm, :bus)], base_name="vvm",
        lower_bound = ref(pm, :bus, i, "vmin")^2,
        upper_bound = ref(pm, :bus, i, "vmax")^2,
        start = 1.0
    )

    vm_delta = var(pm)[:vm_delta] = @variable(pm.model,
        [i in ids(pm, :bus)], base_name="vm_delta",
        start = 0.1
    )

    pg_delta = var(pm)[:pg_delta] = @variable(pm.model,
        [i in ids(pm, :gen)], base_name="pg_delta",
        start = 0.0
    )

    sm_slack = var(pm)[:sm_slack_sparse] = @variable(pm.model,
        [l in ref(pm, :active_rates)], base_name="sm_slack_sparse",
        lower_bound = 0.0,
        start = 0.0
    )

    vr = var(pm, :vr)
    vi = var(pm, :vi)

    PowerModels.constraint_model_voltage(pm)

    for (i,bus) in ref(pm, :bus)
        vm_midpoint = (bus["vmax"] + bus["vmin"])/2.0
        vm_target = min(vm_midpoint + 0.04, bus["vmax"])
        #vm_target = bus["vm"]

        @constraint(pm.model, vr[i]^2 + vi[i]^2 == vvm[i])
        @constraint(pm.model, vvm[i] == vm_target^2 + vm_delta[i])
    end

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end
    #Memento.info(LOGGER, "misc constraints time: $(time() - start_time)")


    start_time = time()
    p = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    q = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    for (i,branch) in ref(pm, :branch)
        #PowerModels.constraint_ohms_yt_from(pm, i)
        #PowerModels.constraint_ohms_yt_to(pm, i)

        f_bus_id = branch["f_bus"]
        t_bus_id = branch["t_bus"]
        f_idx = (i, f_bus_id, t_bus_id)
        t_idx = (i, t_bus_id, f_bus_id)

        f_bus = ref(pm, :bus, f_bus_id)
        t_bus = ref(pm, :bus, t_bus_id)

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        #g = branch["g"]
        #b = branch["b"]
        #tr = branch["tr"]
        #ti = branch["ti"]

        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        #p_fr  = var(pm, :p, f_idx)
        #q_fr  = var(pm, :q, f_idx)
        #p_to  = var(pm, :p, t_idx)
        #q_to  = var(pm, :q, t_idx)
        vr_fr = vr[f_bus_id] #var(pm, :vr, f_bus_id)
        vr_to = vr[t_bus_id] #var(pm, :vr, t_bus_id)
        vi_fr = vi[f_bus_id] #var(pm, :vi, f_bus_id)
        vi_to = vi[t_bus_id] #var(pm, :vi, t_bus_id)


        if branch["transformer"]
            p[f_idx] = (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        else
            p[f_idx] = (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end
        p[t_idx] = (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))

        if haskey(branch, "rate_a")
            #@constraint(pm.model,  p[f_idx] <= branch["rate_a"] + sm_slack[i])
            #@constraint(pm.model, -p[f_idx] <= branch["rate_a"] + sm_slack[i])
            #@constraint(pm.model,  p[t_idx] <= branch["rate_a"] + sm_slack[i])
            #@constraint(pm.model, -p[t_idx] <= branch["rate_a"] + sm_slack[i])

            p_fr = @variable(pm.model, base_name="p_fr", start = 0.0)
            q_fr = @variable(pm.model, base_name="q_fr", start = 0.0)
            p_to = @variable(pm.model, base_name="p_to", start = 0.0)
            q_to = @variable(pm.model, base_name="q_to", start = 0.0)

            @constraint(pm.model,  p[f_idx] == p_fr)
            @constraint(pm.model,  q[f_idx] == q_fr)
            @constraint(pm.model,  p[t_idx] == p_to)
            @constraint(pm.model,  q[t_idx] == q_to)

            @constraint(pm.model, p_fr^2 + q_fr^2 <= (branch["rate_a"] + sm_slack[i])^2)
            @constraint(pm.model, p_to^2 + q_to^2 <= (branch["rate_a"] + sm_slack[i])^2)
        end
    end
    #Memento.info(LOGGER, "flow expr time: $(time() - start_time)")


    start_time = time()
    pg = var(pm, :pg)
    qg = var(pm, :qg)
    for (i,gen) in ref(pm, :gen)
        @constraint(pm.model, pg[i] == gen["pg"] + pg_delta[i])
    end
    #Memento.info(LOGGER, "gen expr time: $(time() - start_time)")


    bsh = var(pm, :bsh)
    for (i,bus) in ref(pm, :bus)
        #PowerModels.constraint_power_balance_shunt(pm, i)

        bus_arcs = ref(pm, :bus_arcs, i)
        bus_arcs_dc = ref(pm, :bus_arcs_dc, i)
        bus_gens = ref(pm, :bus_gens, i)
        bus_loads = ref(pm, :bus_loads, i)
        #bus_shunts = ref(pm, :bus_shunts, i)
        bus_shunts_const = ref(pm, :bus_shunts_const, i)
        bus_shunts_var = ref(pm, :bus_shunts_var, i)

        bus_pd = Dict(k => ref(pm, :load, k, "pd") for k in bus_loads)
        bus_qd = Dict(k => ref(pm, :load, k, "qd") for k in bus_loads)

        bus_gs_const = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts_const)
        bus_bs_const = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts_const)

        #p = var(pm, :p)
        #q = var(pm, :q)
        #pg = var(pm, :pg)
        #qg = var(pm, :qg)

        #@constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr[i]^2 + vi[i]^2))
        #@constraint(pm.model, q_delta[i] + sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr[i]^2 + vi[i]^2))
        @constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vvm[i])
        @constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vvm[i] + sum(bsh[s]*vvm[i] for s in bus_shunts_var))
        #@constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*(vi[i]^2 + vr[i]^2) + sum(bsh[s]*1.0 for s in bus_shunts_var))
    end
    #Memento.info(LOGGER, "power balance constraint time: $(time() - start_time)")

    @objective(pm.model, Min,
        sum( 1e7*vm_delta[i]^2 for (i,bus) in ref(pm, :bus)) +
        sum( 5e5*sm_slack[i] for i in ref(pm, :active_rates)) +
        sum( 1e5*pg_delta[i]^2 for (i,gen) in ref(pm, :gen))
    )
end
