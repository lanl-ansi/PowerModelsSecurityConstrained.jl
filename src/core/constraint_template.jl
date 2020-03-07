""
function constraint_thermal_limit_from_soft(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        constraint_thermal_limit_from_soft(pm, nw, f_idx, branch["rate_a"])
    end
end


""
function constraint_thermal_limit_to_soft(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        constraint_thermal_limit_to_soft(pm, nw, t_idx, branch["rate_a"])
    end
end


""
function constraint_gen_active_deviation(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    gen = ref(pm, nw, :gen, i)

    constraint_gen_active_deviation(pm, nw, i, gen["pg"])
end



""
function constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    cut = ref(pm, :branch_flow_cuts, i)
    branch = ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end

""
function constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    cut = ref(pm, :branch_flow_cuts, i)
    branch = ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end

# revisit after issue #14 is closed
# ""
# function constraint_gen_contingency_ptdf_thermal_limit_from_soft(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
#     cut = ref(pm, :gen_flow_cuts, i)
#     branch = ref(pm, nw, :branch, cut.branch_id)
#     gen = ref(pm, :gen, cut.gen_id)
#     gen_bus = ref(pm, :bus, gen["gen_bus"])
#     gen_set = Set(i for i in ref(pm, :area_gens)[gen_bus["area"]] if i != cut.gen_id)

#     gen_alpha = Dict(i => ref(pm, :gen, i)["alpha"] for i in gen_set)
#     gen_bus = Dict(i => ref(pm, :gen, i)["gen_bus"] for i in gen_set)

#     if haskey(branch, "rate_c")
#         constraint_gen_contingency_ptdf_thermal_limit_from_soft(pm, nw, i, cut.bus_injection, branch["rate_c"], gen_set, gen_alpha, gen_bus)
#     end
# end


""
function constraint_gen_contingency_ptdf_bus_thermal_limit_from_soft(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    cut = ref(pm, :gen_flow_cuts, i)
    branch = ref(pm, nw, :branch, cut.branch_id)
    gen = ref(pm, :gen, cut.gen_id)

    if haskey(branch, "rate_c")
        constraint_gen_contingency_ptdf_thermal_limit_from_soft(pm, nw, i, cut.bus_injection, branch["rate_c"], gen_set, gen_alpha, gen_bus)
    end
end


