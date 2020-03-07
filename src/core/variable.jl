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


function variable_vvm_delta(pm::AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    vvm_delta = var(pm, nw)[:vvm_delta] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_vvm_delta",
        start = comp_start_value(ref(pm, nw, :bus, i), "vvm_delta_start", 0.1)
    )

    report && sol_component_value(pm, nw, :bus, :vvm_delta, ids(pm, nw, :bus), vvm_delta)
end


function variable_pg_delta(pm::AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    pg_delta = var(pm, nw)[:pg_delta] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg_delta",
        start = comp_start_value(ref(pm, nw, :bus, i), "pg_delta_start", 0.1)
    )

    report && sol_component_value(pm, nw, :gen, :pg_delta, ids(pm, nw, :gen), pg_delta)
end


function variable_branch_contigency_flow_violation(pm::AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    branch_cut_vio = var(pm, nw)[:branch_cut_vio] = JuMP.@variable(pm.model,
        [i in 1:length(ref(pm, :branch_flow_cuts))], base_name="$(nw)_cont_branch_flow_vio",
        #start = comp_start_value(ref(pm, nw, :bus, i), "cont_branch_vio_start")
    )

    if bounded
        for i in 1:length(ref(pm, :branch_flow_cuts))
            JuMP.set_lower_bound(branch_cut_vio[i], 0.0)
        end
    end

    #report && sol_component_value(pm, nw, :gen, :pg_delta, ids(pm, nw, :gen), pg_delta)
end


function variable_gen_contigency_flow_violation(pm::AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    gen_cut_vio = var(pm, nw)[:gen_cut_vio] = JuMP.@variable(pm.model,
        [i in 1:length(ref(pm, :gen_flow_cuts))], base_name="$(nw)_cont_gen_flow_vio",
        #start = comp_start_value(ref(pm, nw, :bus, i), "cont_branch_vio_start")
    )

    if bounded
        for i in 1:length(ref(pm, :gen_flow_cuts))
            JuMP.set_lower_bound(gen_cut_vio[i], 0.0)
        end
    end

    #report && sol_component_value(pm, nw, :gen, :pg_delta, ids(pm, nw, :gen), pg_delta)
end


function variable_gen_contigency_capacity_violation(pm::AbstractPowerModel; nw::Int=pm.cnw, bounded::Bool=true, report::Bool=true)
    gen_cap_vio = var(pm, nw)[:gen_cap_vio] = JuMP.@variable(pm.model,
        [i in 1:length(ref(pm, :gen_contingencies))], base_name="$(nw)_cont_gen_cap_vio",
        #start = comp_start_value(ref(pm, nw, :bus, i), "cont_branch_vio_start")
    )

    if bounded
        for i in 1:length(ref(pm, :gen_contingencies))
            JuMP.set_lower_bound(gen_cap_vio[i], 0.0)
        end
    end

    #report && sol_component_value(pm, nw, :gen, :pg_delta, ids(pm, nw, :gen), pg_delta)
end

