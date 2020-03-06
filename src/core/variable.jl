function variable_branch_flow_slack(pm::AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    sm_slack = var(pm, nw)[:sm_slack] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch_sm_active)], base_name="$(nw)_sm_slack",
        start = comp_start_value(ref(pm, nw, :branch, l), "sm_slack_start")
    )

    if bounded
        for (l,branch) in ref(pm, nw, :branch_sm_active)
            JuMP.set_lower_bound(sm_slack[l], 0.0)
        end
    end

    report && sol_component_value(pm, nw, :branch, :sm_slack, ids(pm, nw, :branch_sm_active), sm_slack)
end
