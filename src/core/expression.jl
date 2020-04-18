""
function expression_bus_generation(pm::_PM.AbstractPowerModel, n::Int, i::Int, bus_gens)
    pg = get(var(pm, n), :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg = get(var(pm, n), :qg, Dict()); _PM._check_var_keys(pg, bus_gens, "reactive power", "generator")

    pg_total = 0.0
    qg_total = 0.0
    if length(bus_gens) > 0
        pg_total = sum(pg[g] for g in bus_gens)
        qg_total = sum(qg[g] for g in bus_gens)
    end

    var(pm, n, :bus_pg)[i] = pg_total
    var(pm, n, :bus_qg)[i] = qg_total
end