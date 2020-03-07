""
function constraint_thermal_limit_from_soft(pm::AbstractActivePowerModel, n::Int, f_idx, rate_a)
    l,i,j = f_idx
    p_fr = var(pm, n, :p, f_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_fr <= rate_a + sm_slack)
end

""
function constraint_thermal_limit_to_soft(pm::AbstractActivePowerModel, n::Int, t_idx, rate_a)
    l,i,j = t_idx
    p_to = var(pm, n, :p, t_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_to <= rate_a + sm_slack)
end


""
function expression_bus_generation(pm::AbstractActivePowerModel, n::Int, i::Int, bus_gens)
    pg = get(var(pm, n), :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")

    pg_total = 0.0
    if length(bus_gens) > 0
        pg_total = sum(pg[g] for g in bus_gens)
    end

    var(pm, n, :bus_pg)[i] = pg_total
end


""
function expression_bus_withdrawal(pm::AbstractActivePowerModel, n::Int, i::Int, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    ps = get(var(pm, n), :ps, Dict()); PowerModels._check_var_keys(ps, bus_storage, "active power", "storage")

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








