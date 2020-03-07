""
function constraint_thermal_limit_from_soft(pm::AbstractPowerModel, n::Int, f_idx, rate_a)
    l,i,j = f_idx
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (rate_a + sm_slack)^2)
end


""
function constraint_thermal_limit_to_soft(pm::AbstractPowerModel, n::Int, t_idx, rate_a)
    l,i,j = t_idx
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    sm_slack = var(pm, :sm_slack, l)

    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (rate_a + sm_slack)^2)
end


""
function constraint_gen_active_deviation(pm::AbstractPowerModel, n::Int, i, pg)
    pg_var = var(pm, n, :pg, i)
    pg_delta = var(pm, n, :pg_delta, i)

    JuMP.@constraint(pm.model, pg_var == pg + pg_delta)
end