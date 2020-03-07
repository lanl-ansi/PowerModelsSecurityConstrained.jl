
""
function constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm::AbstractDCPModel, n::Int, i::Int, cut_map, rate)
    bus_injection = var(pm, :bus_pg)
    bus_withdrawal = var(pm, :bus_wdp)

    @constraint(pm.model, sum(weight*(bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut_map) <= rate + var(pm, :branch_cut_vio, i))
end


""
function constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm::AbstractDCPModel, n::Int, i::Int, cut_map, rate)
    bus_injection = var(pm, :bus_pg)
    bus_withdrawal = var(pm, :bus_wdp)

    @constraint(pm.model, -sum(weight*(bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut_map) <= rate + var(pm, :branch_cut_vio, i))
end


# revisit after issue #14 is closed
# ""
# function constraint_gen_contingency_ptdf_thermal_limit_from_soft(pm::AbstractDCPModel, n::Int, i::Int, cut_map, rate, gen_set, gen_alpha, gen_bus)
# end

