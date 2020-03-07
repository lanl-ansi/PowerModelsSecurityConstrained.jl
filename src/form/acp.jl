

""
function constraint_power_balance_shunt_dispatch(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vm = var(pm, n, :vm, i)
    p = var(pm, n, :p)
    q = var(pm, n, :q)
    pg = var(pm, n, :pg)
    qg = var(pm, n, :qg)
    bs = var(pm, n, :bs)
    p_dc = var(pm, n, :p_dc)
    q_dc = var(pm, n, :q_dc)

    cstr_p = JuMP.@NLconstraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm^2)
    cstr_q = JuMP.@NLconstraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm^2 + sum(bs[s]*vm^2 for s in bus_shunts_var))

    if report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

""
function constraint_power_balance_shunt_dispatch_soft(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vm = var(pm, n, :vm, i)
    p_delta_abs = var(pm, n, :p_delta_abs, i)
    q_delta_abs = var(pm, n, :q_delta_abs, i)

    p = var(pm, n, :p)
    q = var(pm, n, :q)
    pg = var(pm, n, :pg)
    qg = var(pm, n, :qg)
    bs = var(pm, n, :bs)
    p_dc = var(pm, n, :p_dc)
    q_dc = var(pm, n, :q_dc)

    p_delta = @NLexpression(pm.model, - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm^2)
    q_delta = @NLexpression(pm.model, - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm^2 + sum(bs[s]*vm^2 for s in bus_shunts_var))

    @NLconstraint(pm.model,  p_delta_abs >= p_delta)
    @NLconstraint(pm.model, -p_delta_abs <= p_delta)

    @NLconstraint(pm.model,  q_delta_abs >= q_delta)
    @NLconstraint(pm.model, -q_delta_abs <= q_delta)
end


function constraint_ohms_yt_from_goc(pm::AbstractACPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n,  :p, f_idx)
    q_fr  = var(pm, n,  :q, f_idx)
    vm_fr = var(pm, n, :vm, f_bus)
    vm_to = var(pm, n, :vm, t_bus)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_fr ==  (g/tm^2+g_fr)*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b/tm^2+b_fr)*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
end


""
function expression_bus_withdrawal(pm::AbstractACPModel, n::Int, i::Int, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm = var(pm, n, :vm, i)
    ps = get(var(pm, n), :ps, Dict()); PowerModels._check_var_keys(ps, bus_storage, "active power", "storage")
    qs = get(var(pm, n), :qs, Dict()); PowerModels._check_var_keys(ps, bus_storage, "reactive power", "storage")

    ps_total = 0.0
    qs_total = 0.0
    if length(bus_storage) > 0
        ps_total = sum(ps[s] for s in bus_storage)
        qs_total = sum(qs[s] for s in bus_storage)
    end

    pd_total = 0.0
    qd_total = 0.0
    if length(bus_pd) > 0
        pd_total = sum(pd for pd in values(bus_pd))
    end
        if length(bus_qd) > 0
        qd_total = sum(qd for qd in values(bus_qd))
    end

    gs_total = 0.0
    bs_total = 0.0
    if length(bus_gs) > 0
        gs_total = sum(gs for gs in values(bus_gs))*vm^2
    end
    if length(bus_bs) > 0
        bs_total = sum(bs for bs in values(bus_bs))*vm^2
    end

    var(pm, n, :bus_wdp)[i] = ps_total + pd_total + gs_total
    var(pm, n, :bus_wdq)[i] = qs_total + qd_total + bs_total
end

