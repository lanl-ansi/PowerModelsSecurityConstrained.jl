
""
function variable_c1_shunt_admittance_imaginary(pm::_PM.AbstractActivePowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
# do nothing
end


""
function constraint_c1_power_balance_shunt_dispatch(pm::_PM.AbstractActivePowerModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    p    = get(var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _PM._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")

    cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*1.0)

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
    end
end


""
function constraint_c1_thermal_limit_from_soft(pm::_PM.AbstractActivePowerModel, n::Int, f_idx, rate_a)
    l,i,j = f_idx
    p_fr = var(pm, n, :p, f_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_fr <= rate_a + sm_slack)
end

""
function constraint_c1_thermal_limit_to_soft(pm::_PM.AbstractActivePowerModel, n::Int, t_idx, rate_a)
    l,i,j = t_idx
    p_to = var(pm, n, :p, t_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_to <= rate_a + sm_slack)
end


""
function expression_c1_bus_generation(pm::_PM.AbstractActivePowerModel, n::Int, i::Int, bus_gens)
    pg = get(var(pm, n), :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")

    pg_total = 0.0
    if length(bus_gens) > 0
        pg_total = sum(pg[g] for g in bus_gens)
    end

    var(pm, n, :bus_pg)[i] = pg_total
end


""
function expression_c1_bus_withdrawal(pm::_PM.AbstractActivePowerModel, n::Int, i::Int, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    ps = get(var(pm, n), :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")

    ps_total = 0.0
    if length(bus_storage) > 0
        ps_total = sum(ps[s] for s in bus_storage)
    end

    pd_total = 0.0
    if length(bus_pd) > 0
        pd_total = sum(pd for pd in values(bus_pd))
    end

    gs_total = 0.0
    if length(bus_gs) > 0
        gs_total = sum(gs for gs in values(bus_gs))*1.0^2
    end

    var(pm, n, :bus_wdp)[i] = ps_total + pd_total + gs_total
end





""
function constraint_c2_power_balance(pm::_PM.AbstractActivePowerModel, n::Int, i::Int, bus_arcs, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
    p    = get(var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")

    z_demand = get(var(pm, n), :z_demand, Dict()); _PM._check_var_keys(z_demand, keys(bus_pd), "power factor", "load")

    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(pd*z_demand[i] for (i,pd) in bus_pd)
        - sum(gs*1.0 for (i,gs) in bus_gs)*1.0^2
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
    end
end

""
function constraint_c2_power_balance_soft_lin(pm::_PM.AbstractActivePowerModel, n::Int, i::Int, bus_arcs, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
    p_delta_abs = var(pm, n, :p_delta_abs, i)
    p    = get(var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")

    z_demand = get(var(pm, n), :z_demand, Dict()); _PM._check_var_keys(z_demand, keys(bus_pd), "power factor", "load")

    JuMP.@constraint(pm.model,
        p_delta_abs >=
        sum(pg[g] for g in bus_gens)
        - sum(pd*z_demand[i] for (i,pd) in bus_pd)
        - sum(gs*1.0 for (i,gs) in bus_gs)*1.0^2
        - sum(p[a] for a in bus_arcs)
    )
    JuMP.@constraint(pm.model,
        p_delta_abs >= -(
        sum(pg[g] for g in bus_gens)
        - sum(pd*z_demand[i] for (i,pd) in bus_pd)
        - sum(gs*1.0 for (i,gs) in bus_gs)*1.0^2
        - sum(p[a] for a in bus_arcs))
    )
end

""
function constraint_c2_current_limit_from_soft(pm::_PM.AbstractActivePowerModel, n::Int, f_idx, rate_a)
    l,i,j = f_idx
    p_fr = var(pm, n, :p, f_idx)
    sm_slack = var(pm, n, :sm_slack, l)

    JuMP.@constraint(pm.model, p_fr <= rate_a*(1 + sm_slack))
end

""
function constraint_c2_current_limit_to_soft(pm::_PM.AbstractActivePowerModel, n::Int, t_idx, rate_a)
    l,i,j = t_idx
    p_to = var(pm, n, :p, t_idx)
    sm_slack = var(pm, n, :sm_slack, l)

    JuMP.@constraint(pm.model, p_to <= rate_a*(1 + sm_slack))
end


""
function constraint_c2_thermal_limit_from_soft(pm::_PM.AbstractActivePowerModel, n::Int, f_idx, rate_a)
    l,i,j = f_idx
    p_fr = var(pm, n, :p, f_idx)
    sm_slack = var(pm, n, :sm_slack, l)

    JuMP.@constraint(pm.model, p_fr <= rate_a*(1 + sm_slack))
end

""
function constraint_c2_thermal_limit_to_soft(pm::_PM.AbstractActivePowerModel, n::Int, t_idx, rate_a)
    l,i,j = t_idx
    p_to = var(pm, n, :p, t_idx)
    sm_slack = var(pm, n, :sm_slack, l)

    JuMP.@constraint(pm.model, p_to <= rate_a*(1 + sm_slack))
end




