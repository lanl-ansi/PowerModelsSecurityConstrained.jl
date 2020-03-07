
""
function constraint_power_balance_shunt_dispatch(pm::AbstractACRModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vi = var(pm, n, :vi, i)
    vr = var(pm, n, :vr, i)
    p = var(pm, n, :p)
    q = var(pm, n, :q)
    pg = var(pm, n, :pg)
    qg = var(pm, n, :qg)
    bs = var(pm, n, :bs)
    p_dc = var(pm, n, :p_dc)
    q_dc = var(pm, n, :q_dc)

    # possibly can save 2x in function eval, but no the dominant runtime at this moment
    #vm_sqr = @variable(pm.model, start=1.0, base_name="$(0)_vm_sqr_$(i)")

    #JuMP.@constraint(pm.model, vm_sqr == vi^2 + vr^2)
    #cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm_sqr)
    #cstr_q = JuMP.@constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm_sqr + sum(bs[s]*vm_sqr for s in bus_shunts_var))

    cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*(vi^2 + vr^2))
    cstr_q = JuMP.@NLconstraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*(vi^2 + vr^2) + sum(bs[s]*(vi^2 + vr^2) for s in bus_shunts_var))

    if report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


function constraint_ohms_yt_from_goc(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    vr_fr = var(pm, n, :vr, f_bus)
    vr_to = var(pm, n, :vr, t_bus)
    vi_fr = var(pm, n, :vi, f_bus)
    vi_to = var(pm, n, :vi, t_bus)

    JuMP.@constraint(pm.model, p_fr ==  (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
    JuMP.@constraint(pm.model, q_fr == -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
end
