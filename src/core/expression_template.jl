"""
defines power from generators at each bus, varies in contingencies
"""
function expression_c1_bus_generation(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(var(pm, nw), :bus_pg)
        var(pm, nw)[:bus_pg] = Dict{Int,Any}()
    end
    if !haskey(var(pm, nw), :bus_qg)
        var(pm, nw)[:bus_qg] = Dict{Int,Any}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_gens = ref(pm, nw, :bus_gens, i)

    expression_c1_bus_generation(pm, nw, i, bus_gens)
end


"""
defines power from non-generator components at each bus, static in contingencies
"""
function expression_c1_bus_withdrawal(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(var(pm, nw), :bus_wdp)
        var(pm, nw)[:bus_wdp] = Dict{Int,Any}()
    end
    if !haskey(var(pm, nw), :bus_wdq)
        var(pm, nw)[:bus_wdq] = Dict{Int,Any}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    expression_c1_bus_withdrawal(pm, nw, i, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end



""
function expression_c1_branch_power_ohms_yt_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(var(pm, nw), :p)
        var(pm, nw)[:p] = Dict{Tuple{Int,Int,Int},Any}()
    end
    if !haskey(var(pm, nw), :q)
        var(pm, nw)[:q] = Dict{Tuple{Int,Int,Int},Any}()
    end

    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    if branch["transformer"]
        expression_c1_branch_power_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    else
        _PM.expression_branch_power_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    end
end



