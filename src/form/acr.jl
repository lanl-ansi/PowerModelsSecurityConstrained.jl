
""
function constraint_power_balance_shunt_dispatch(pm::_PM.AbstractACRModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vi = var(pm, n, :vi, i)
    vr = var(pm, n, :vr, i)
    p    = get(var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, n),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); _PM._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); _PM._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    bs = get(var(pm, n), :bs, Dict()); _PM._check_var_keys(bs, bus_shunts_var, "reactive power", "shunt")


    # possibly can save 2x in function eval, but no the dominant runtime at this moment
    #vm_sqr = @variable(pm.model, start=1.0, base_name="$(0)_vm_sqr_$(i)")

    #JuMP.@constraint(pm.model, vm_sqr == vi^2 + vr^2)
    #cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm_sqr)
    #cstr_q = JuMP.@constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm_sqr + sum(bs[s]*vm_sqr for s in bus_shunts_var))

    cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*(vi^2 + vr^2))
    cstr_q = JuMP.@NLconstraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*(vi^2 + vr^2) + sum(bs[s]*(vi^2 + vr^2) for s in bus_shunts_var))

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


function constraint_ohms_yt_from_goc(pm::_PM.AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    JuMP.@constraint(pm.model, p_fr ==  (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
    JuMP.@constraint(pm.model, q_fr == -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
end




""
function expression_branch_power_ohms_yt_from_goc(pm::_PM.AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    var(pm, n, :p)[f_idx] =  (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
    var(pm, n, :q)[f_idx] = -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
end

