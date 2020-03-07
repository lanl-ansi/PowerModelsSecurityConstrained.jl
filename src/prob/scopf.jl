

"""
An SCOPF formulation conforming to the ARPA-e GOC Challenge 1 specification.
A DC power flow approximation is used. Power balance and line flow constraints
are strictly enforced in the first stage.  Contingency branch flow constraints
are enforced by PTDF cuts and penalized based on a conservative linear
approximation of the formulation's specification.

This formulation is used in conjunction with the contingency filters that
generate PTDF cuts.
"""
function run_scopf_cuts_soft(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, build_scopf_cuts_soft; kwargs...)
end

""
function build_scopf_cuts_soft(pm::AbstractPowerModel)
    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm)
    PowerModels.variable_dcline_flow(pm)

    variable_branch_contigency_flow_violation(pm)
    variable_gen_contigency_flow_violation(pm)
    variable_gen_contigency_capacity_violation(pm)

    for i in ids(pm, :bus)
        expression_bus_generation(pm, i)
        expression_bus_withdrawal(pm, i)
    end


    PowerModels.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        PowerModels.constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        PowerModels.constraint_ohms_yt_from(pm, i)
        PowerModels.constraint_ohms_yt_to(pm, i)

        PowerModels.constraint_voltage_angle_difference(pm, i)

        PowerModels.constraint_thermal_limit_from(pm, i)
        PowerModels.constraint_thermal_limit_to(pm, i)
    end


    for (i,cut) in enumerate(ref(pm, :branch_flow_cuts))
        constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm, i)
        constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, i)
    end

    bus_withdrawal = var(pm, :bus_wdp)

    for (i,cut) in enumerate(ref(pm, :gen_flow_cuts))
        branch = ref(pm, :branch, cut.branch_id)
        gen = ref(pm, :gen, cut.gen_id)
        gen_bus = ref(pm, :bus, gen["gen_bus"])
        gen_set = ref(pm, :area_gens)[gen_bus["area"]]
        alpha_total = sum(gen["alpha"] for (i,gen) in ref(pm, :gen) if gen["index"] != cut.gen_id && i in gen_set)

        cont_bus_injection = Dict{Int,Any}()
        for (i, bus) in ref(pm, :bus)
            inj = 0.0
            for g in ref(pm, :bus_gens, i)
                if g != cut.gen_id
                    if g in gen_set
                        inj += var(pm, :pg, g) + gen["alpha"]*var(pm, :pg, cut.gen_id)/alpha_total
                    else
                        inj += var(pm, :pg, g)
                    end
                end
            end
            cont_bus_injection[i] = inj
        end

        #rate = branch["rate_a"]
        rate = branch["rate_c"]
        @constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + var(pm, :gen_cont_flow_vio, i))
        @constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + var(pm, :gen_cont_flow_vio, i))
    end

    for (i,gen_cont) in enumerate(ref(pm, :gen_contingencies))
        #println(gen_cont)
        gen = ref(pm, :gen, gen_cont.idx)
        gen_bus = ref(pm, :bus, gen["gen_bus"])
        gen_set = ref(pm, :area_gens)[gen_bus["area"]]
        response_gens = Dict(g => ref(pm, :gen, g) for g in gen_set if g != gen_cont.idx)

        # factor of 1.2 accounts for losses in a DC model
        #@constraint(pm.model, sum(gen["pmax"] - var(pm, :pg, g) for (g,gen) in response_gens) >= 1.2*var(pm, :pg, gen_cont.idx))
        @constraint(pm.model, var(pm, :gen_cont_cap_vio, i) + sum(gen["pmax"] - var(pm, :pg, g) for (g,gen) in response_gens) >= var(pm, :pg, gen_cont.idx))
        #@constraint(pm.model, sum(gen["pmin"] - var(pm, :pg, g) for (g,gen) in response_gens) <= var(pm, :pg, gen_cont.idx))
    end

    ##### Setup Objective #####
    objective_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = var(pm, pm.cnw, :pg_cost)
    branch_cont_flow_vio = var(pm, pm.cnw, :branch_cont_flow_vio)
    gen_cont_flow_vio = var(pm, pm.cnw, :gen_cont_flow_vio)
    gen_cont_cap_vio = var(pm, pm.cnw, :gen_cont_cap_vio)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) ) +
        sum( 5e5*branch_cont_flow_vio[i] for i in 1:length(ref(pm, :branch_flow_cuts)) ) +
        sum( 5e5*gen_cont_flow_vio[i] for i in 1:length(ref(pm, :gen_flow_cuts)) ) + 
        sum( 5e5*gen_cont_cap_vio[i] for i in 1:length(ref(pm, :gen_contingencies)) )
    )
end



"a variant of `run_scopf_cuts_dc_soft` with a different generator response function"
function run_scopf_cuts_soft_2(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, build_scopf_cuts_soft_2; kwargs...)
end

""
function build_scopf_cuts_soft_2(pm::AbstractPowerModel)
    PowerModels.variable_voltage(pm)
    PowerModels.variable_generation(pm)
    PowerModels.variable_branch_flow(pm)
    PowerModels.variable_dcline_flow(pm)

    variable_branch_contigency_flow_violation(pm)
    variable_gen_contigency_flow_violation(pm)

    for i in ids(pm, :bus)
        expression_bus_generation(pm, i)
        expression_bus_withdrawal(pm, i)
    end


    PowerModels.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        PowerModels.constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        PowerModels.constraint_ohms_yt_from(pm, i)
        PowerModels.constraint_ohms_yt_to(pm, i)

        PowerModels.constraint_voltage_angle_difference(pm, i)

        PowerModels.constraint_thermal_limit_from(pm, i)
        PowerModels.constraint_thermal_limit_to(pm, i)
    end

    for (i,cut) in enumerate(ref(pm, :branch_flow_cuts))
        constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm, i)
        constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, i)
    end


    bus_withdrawal = var(pm, :bus_wdp)

    for (i,cut) in enumerate(ref(pm, :gen_flow_cuts))
        branch = ref(pm, :branch, cut.branch_id)

        bus_num = length(ref(pm, :bus))
        cont_bus_injection = Dict{Int,Any}()
        for (b,bus) in ref(pm, :bus)
            inj = 0.0
            for g in ref(pm, :bus_gens, b)
                if g != cut.gen_id
                    inj += var(pm, :pg, g)
                end
            end
            cont_bus_injection[b] = inj + var(pm, :pg, cut.gen_id)/bus_num
        end

        #rate = branch["rate_a"]
        rate = branch["rate_c"]
        @constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + var(pm, :gen_cont_flow_vio, i))
        @constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + var(pm, :gen_cont_flow_vio, i))
        #@constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate)
        #@constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate)
    end

    ##### Setup Objective #####
    objective_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = var(pm, pm.cnw, :pg_cost)
    branch_cont_flow_vio = var(pm, pm.cnw, :branch_cont_flow_vio)
    gen_cont_flow_vio = var(pm, pm.cnw, :gen_cont_flow_vio)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) ) +
        sum( 5e5*branch_cont_flow_vio[i] for i in 1:length(ref(pm, :branch_flow_cuts)) ) +
        sum( 5e5*gen_cont_flow_vio[i] for i in 1:length(ref(pm, :gen_flow_cuts)) )
    )
end
