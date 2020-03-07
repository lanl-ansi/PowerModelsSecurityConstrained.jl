"""
defines power from generators at each bus, varies in contingencies
"""
function expression_bus_generation(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    if !haskey(var(pm, nw), :bus_pg)
        var(pm, nw)[:bus_pg] = Dict{Int,Any}()
    end
    if !haskey(var(pm, nw), :bus_qg)
        var(pm, nw)[:bus_qg] = Dict{Int,Any}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_gens = ref(pm, nw, :bus_gens, i)

    expression_bus_generation(pm, nw, i, bus_gens)
end


"""
defines power from non-generator components at each bus, static in contingencies
"""
function expression_bus_withdrawal(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
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

    expression_bus_withdrawal(pm, nw, i, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end
