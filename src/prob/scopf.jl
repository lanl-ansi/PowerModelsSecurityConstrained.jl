"""
An academic SCOPF formulation inspired by the ARPA-e GOC Challenge 1 specification.
Power balance and line flow constraints are strictly enforced in the first
stage and contingency stages.

This formulation is best used in conjunction with the contingency filters that
find violated contingencies.
"""
function run_c1_scopf(file, model_constructor, solver; kwargs...)
    return _PM.solve_model(file, model_constructor, solver, build_c1_scopf; multinetwork=true, kwargs...)
end

# enables support for v[1], required for objective_variable_pg_cost when pg is an expression
Base.getindex(v::JuMP.GenericAffExpr, i::Int64) = v

""
function build_c1_scopf(pm::_PM.AbstractPowerModel)
    # base-case network id is 0

    _PM.variable_bus_voltage(pm, nw=0)
    _PM.variable_gen_power(pm, nw=0)
    _PM.variable_branch_power(pm, nw=0)

    _PM.constraint_model_voltage(pm, nw=0)

    for i in ids(pm, :ref_buses, nw=0)
        _PM.constraint_theta_ref(pm, i, nw=0)
    end

    for i in ids(pm, :bus, nw=0)
        _PM.constraint_power_balance(pm, i, nw=0)
    end

    for i in ids(pm, :branch, nw=0)
        _PM.constraint_ohms_yt_from(pm, i, nw=0)
        _PM.constraint_ohms_yt_to(pm, i, nw=0)

        _PM.constraint_voltage_angle_difference(pm, i, nw=0)

        _PM.constraint_thermal_limit_from(pm, i, nw=0)
        _PM.constraint_thermal_limit_to(pm, i, nw=0)
    end


    contigency_ids = [id for id in nw_ids(pm) if id != 0]
    for nw in contigency_ids
        _PM.variable_bus_voltage(pm, nw=nw, bounded=false)
        _PM.variable_gen_power(pm, nw=nw, bounded=false)
        _PM.variable_branch_power(pm, nw=nw)

        variable_c1_response_delta(pm, nw=nw)


        _PM.constraint_model_voltage(pm, nw=nw)

        for i in ids(pm, :ref_buses, nw=nw)
            _PM.constraint_theta_ref(pm, i, nw=nw)
        end

        gen_buses = ref(pm, :gen_buses, nw=nw)
        for i in ids(pm, :bus, nw=nw)
            _PM.constraint_power_balance(pm, i, nw=nw)

            # if a bus has active generators, fix the voltage magnitude to the base case
            if i in gen_buses
                constraint_c1_voltage_magnitude_link(pm, i, nw_1=0, nw_2=nw)
            end
        end


        response_gens = ref(pm, :response_gens, nw=nw)
        for (i,gen) in ref(pm, :gen, nw=nw)
            pg_base = var(pm, :pg, i, nw=0)

            # setup the linear response function or fix value to base case
            if i in response_gens
                constraint_c1_gen_power_real_response(pm, i, nw_1=0, nw_2=nw)
            else
                constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)
            end
        end


        for i in ids(pm, :branch, nw=nw)
            _PM.constraint_ohms_yt_from(pm, i, nw=nw)
            _PM.constraint_ohms_yt_to(pm, i, nw=nw)

            _PM.constraint_voltage_angle_difference(pm, i, nw=nw)

            _PM.constraint_thermal_limit_from(pm, i, nw=nw)
            _PM.constraint_thermal_limit_to(pm, i, nw=nw)
        end
    end


    ##### Setup Objective #####
    objective_c1_variable_pg_cost_basecase(pm)

    # explicit network id needed because of conductor-less
    pg_cost = var(pm, 0, :pg_cost)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, 0, :gen) )
    )
end



"""
An academic SCOPF formulation inspired by the ARPA-e GOC Challenge 1 specification.
Power balance and line flow constraints are strictly enforced in the first
stage and contingency stages. Contingency branch flow constraints are enforced
by PTDF cuts using the DC power flow approximation.

This formulation is used in conjunction with the contingency filters that
generate PTDF cuts.
"""
function run_c1_scopf_cuts(file, model_constructor, solver; kwargs...)
    return _PM.solve_model(file, model_constructor, solver, build_c1_scopf_cuts; kwargs...)
end

""
function build_c1_scopf_cuts(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)

    for i in ids(pm, :bus)
        expression_c1_bus_generation(pm, i)
        expression_c1_bus_withdrawal(pm, i)
    end

    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        _PM.constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end


    for (i,cut) in enumerate(ref(pm, :branch_flow_cuts))
        constraint_c1_branch_contingency_ptdf_thermal_limit_from(pm, i)
        constraint_c1_branch_contingency_ptdf_thermal_limit_to(pm, i)
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
        @constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate)
        @constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate)
    end

    ##### Setup Objective #####
    objective_c1_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = var(pm, :pg_cost)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) )
    )
end



"""
An SCOPF formulation conforming to the ARPA-e GOC Challenge 1 specification.
A DC power flow approximation is used. Power balance and line flow constraints
are strictly enforced in the first stage.  Contingency branch flow constraints
are enforced by PTDF cuts and penalized based on a conservative linear
approximation of the formulation's specification.

This formulation is used in conjunction with the contingency filters that
generate PTDF cuts.
"""
function run_c1_scopf_cuts_soft(file, model_constructor, solver; kwargs...)
    return _PM.solve_model(file, model_constructor, solver, build_c1_scopf_cuts_soft; ref_extensions=[ref_c1!], kwargs...)
end

""
function build_c1_scopf_cuts_soft(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)

    variable_c1_shunt_admittance_imaginary(pm)

    variable_c1_branch_contigency_power_violation(pm)
    variable_c1_gen_contigency_power_violation(pm)
    variable_c1_gen_contigency_capacity_violation(pm)

    for i in ids(pm, :bus)
        expression_c1_bus_generation(pm, i)
        expression_c1_bus_withdrawal(pm, i)
    end


    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_c1_power_balance_shunt_dispatch(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end


    for (i,cut) in enumerate(ref(pm, :branch_flow_cuts))
        constraint_c1_branch_contingency_ptdf_thermal_limit_from_soft(pm, i)
        constraint_c1_branch_contingency_ptdf_thermal_limit_to_soft(pm, i)
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
    objective_c1_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = var(pm, :pg_cost)
    branch_cont_flow_vio = var(pm, :branch_cont_flow_vio)
    gen_cont_flow_vio = var(pm, :gen_cont_flow_vio)
    gen_cont_cap_vio = var(pm, :gen_cont_cap_vio)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) ) +
        sum( 5e5*branch_cont_flow_vio[i] for i in 1:length(ref(pm, :branch_flow_cuts)) ) +
        sum( 5e5*gen_cont_flow_vio[i] for i in 1:length(ref(pm, :gen_flow_cuts)) ) + 
        sum( 5e5*gen_cont_cap_vio[i] for i in 1:length(ref(pm, :gen_contingencies)) )
    )
end



"a variant of `run_c1_scopf_cuts_dc_soft` with a different generator response function"
function run_c1_scopf_cuts_soft_bpv(file, model_constructor, solver; kwargs...)
    return _PM.solve_model(file, model_constructor, solver, build_c1_scopf_cuts_soft_bpv; ref_extensions=[ref_c1!], kwargs...)
end

""
function build_c1_scopf_cuts_soft_bpv(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)

    variable_c1_shunt_admittance_imaginary(pm)

    variable_c1_branch_contigency_power_violation(pm)
    variable_c1_gen_contigency_power_violation(pm)

    for i in ids(pm, :bus)
        expression_c1_bus_generation(pm, i)
        expression_c1_bus_withdrawal(pm, i)
    end


    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_c1_power_balance_shunt_dispatch(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    for (i,cut) in enumerate(ref(pm, :branch_flow_cuts))
        constraint_c1_branch_contingency_ptdf_thermal_limit_from_soft(pm, i)
        constraint_c1_branch_contingency_ptdf_thermal_limit_to_soft(pm, i)
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
    objective_c1_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = var(pm, :pg_cost)
    branch_cont_flow_vio = var(pm, :branch_cont_flow_vio)
    gen_cont_flow_vio = var(pm, :gen_cont_flow_vio)

    @objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in ref(pm, :gen) ) +
        sum( 5e5*branch_cont_flow_vio[i] for i in 1:length(ref(pm, :branch_flow_cuts)) ) +
        sum( 5e5*gen_cont_flow_vio[i] for i in 1:length(ref(pm, :gen_flow_cuts)) )
    )
end
