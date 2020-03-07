
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

