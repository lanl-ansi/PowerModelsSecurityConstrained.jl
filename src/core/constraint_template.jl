""
function constraint_power_balance_shunt_dispatch(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_shunts_const = ref(pm, :bus_shunts_const, i)
    bus_shunts_var = ref(pm, :bus_shunts_var, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs_const = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts_const)
    bus_bs_const = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts_const)

    constraint_power_balance_shunt_dispatch(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
end

""
function constraint_power_balance_shunt_dispatch_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_shunts_const = ref(pm, :bus_shunts_const, i)
    bus_shunts_var = ref(pm, :bus_shunts_var, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs_const = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts_const)
    bus_bs_const = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts_const)

    constraint_power_balance_shunt_dispatch_soft(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
end


""
function constraint_ohms_yt_from_goc(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
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
        constraint_ohms_yt_from_goc(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    else
        _PM.constraint_ohms_yt_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    end
end


"links the voltage voltage_magnitude of two networks together"
function constraint_voltage_magnitude_link(pm::_PM.AbstractPowerModel, i::Int; nw_1::Int=nw_id_default, nw_2::Int=nw_id_default)
    constraint_voltage_magnitude_link(pm, nw_1, nw_2, i)
end


"links the voltage voltage_magnitude of two networks together"
function constraint_gen_power_real_link(pm::_PM.AbstractPowerModel, i::Int; nw_1::Int=nw_id_default, nw_2::Int=nw_id_default)
    constraint_gen_power_real_link(pm, nw_1, nw_2, i)
end

"links the generator power of two networks together, with a linear response function"
function constraint_gen_power_real_response(pm::_PM.AbstractPowerModel, i::Int; nw_1::Int=nw_id_default, nw_2::Int=nw_id_default)
    gen = ref(pm, nw_2, :gen, i)
    constraint_gen_power_real_response(pm, nw_1, nw_2, i, gen["alpha"])
end


""
function constraint_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        constraint_thermal_limit_from_soft(pm, nw, f_idx, branch["rate_a"])
    end
end


""
function constraint_thermal_limit_to_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        constraint_thermal_limit_to_soft(pm, nw, t_idx, branch["rate_a"])
    end
end


""
function constraint_gen_power_real_deviation(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)

    constraint_gen_power_real_deviation(pm, nw, i, gen["pg"])
end


""
function constraint_branch_contingency_ptdf_thermal_limit_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cut = ref(pm, :branch_flow_cuts, i)
    branch = ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_from(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end

""
function constraint_branch_contingency_ptdf_thermal_limit_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cut = ref(pm, :branch_flow_cuts, i)
    branch = ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_to(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end


""
function constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cut = ref(pm, :branch_flow_cuts, i)
    branch = ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end

""
function constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cut = ref(pm, :branch_flow_cuts, i)
    branch = ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end

# revisit after issue #14 is closed
# ""
# function constraint_gen_contingency_ptdf_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
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

