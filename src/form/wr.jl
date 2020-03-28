

""
function constraint_power_balance_shunt_dispatch(pm::AbstractWRModels, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    w   = var(pm, n, :w, i)
    p    = get(var(pm, n),    :p, Dict()); PowerModels._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, n),    :q, Dict()); PowerModels._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, n),   :qg, Dict()); PowerModels._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); PowerModels._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); PowerModels._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); PowerModels._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); PowerModels._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); PowerModels._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); PowerModels._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    bs = get(var(pm, n), :bs, Dict()); PowerModels._check_var_keys(bs, bus_shunts_var, "reactive power", "shunt")
    wbs = get(var(pm, n), :wbs, Dict()); PowerModels._check_var_keys(bs, bus_shunts_var, "reactive power", "shunt")

    cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w)
    cstr_q = JuMP.@constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbs[s] for s in bus_shunts_var))

    for s in bus_shunts_var
        _IM.relaxation_product(pm.model, w, bs[s], wbs[s])
    end

    if report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


""
function constraint_power_balance_shunt_dispatch_soft(pm::AbstractWRModels, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    w = var(pm, n, :w, i)
    p_delta_abs = var(pm, n, :p_delta_abs, i)
    q_delta_abs = var(pm, n, :q_delta_abs, i)

    p = var(pm, n, :p)
    q = var(pm, n, :q)
    pg = var(pm, n, :pg)
    qg = var(pm, n, :qg)
    bs = var(pm, n, :bs)
    wbs = var(pm, n, :wbs)
    p_dc = var(pm, n, :p_dc)
    q_dc = var(pm, n, :q_dc)

    #p_delta = - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w
    #q_delta = - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbs[s] for s in bus_shunts_var)

    @constraint(pm.model,  p_delta_abs >= - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w)
    @constraint(pm.model, -p_delta_abs <= - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w)

    @constraint(pm.model,  q_delta_abs >= - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbs[s] for s in bus_shunts_var))
    @constraint(pm.model, -q_delta_abs <= - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbs[s] for s in bus_shunts_var))

    for s in bus_shunts_var
        _IM.relaxation_product(pm.model, w, bs[s], wbs[s])
    end
end


function constraint_ohms_yt_from_goc(pm::AbstractWRModels, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    w_fr = var(pm, n, :w, f_bus)
    wr   = var(pm, n, :wr, (f_bus, t_bus))
    wi   = var(pm, n, :wi, (f_bus, t_bus))

    JuMP.@constraint(pm.model, p_fr ==  (g/tm^2+g_fr)*w_fr + (-g*tr+b*ti)/tm^2*wr + (-b*tr-g*ti)/tm^2*wi )
    JuMP.@constraint(pm.model, q_fr == -(b/tm^2+b_fr)*w_fr - (-b*tr-g*ti)/tm^2*wr + (-g*tr+b*ti)/tm^2*wi )
end


