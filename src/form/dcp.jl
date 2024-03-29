
# revisit after issue #14 is closed
# ""
# function constraint_c1_gen_contingency_ptdf_thermal_limit_from_soft(pm::_PM.AbstractDCPModel, n::Int, i::Int, cut_map, rate, gen_set, gen_alpha, gen_bus)
# end

""
function constraint_c1_voltage_magnitude_link(pm::_PM.AbstractDCPModel, n_1::Int, n_2::Int, i::Int)
    # do nothing voltages are assumed to be the same
end

function constraint_goc_ohms_yt_from(pm::_PM.AbstractDCPModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n,  :p, f_idx)
    va_fr = var(pm, n, :va, f_bus)
    va_to = var(pm, n, :va, t_bus)

    JuMP.@constraint(pm.model, p_fr == -b*(va_fr - va_to))
    # omit reactive constraint
end
