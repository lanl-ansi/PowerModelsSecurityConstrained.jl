""
function constraint_thermal_limit_from_soft(pm::AbstractActivePowerModel, n::Int, f_idx, rate_a)
    l,i,j = f_idx
    p_fr = var(pm, n, :p, f_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_fr <= rate_a + sm_slack)
end

""
function constraint_thermal_limit_to_soft(pm::AbstractActivePowerModel, n::Int, t_idx, rate_a)
    l,i,j = t_idx
    p_to = var(pm, n, :p, t_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_to <= rate_a + sm_slack)
end
