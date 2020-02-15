"given delta, computes the total power response"
function comp_pg_response_total(network, gens::Set{Int}; delta=network["delta"])
    total_pg = 0.0

    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0 
            if gen["index"] in gens
                pg = gen["pg_base"] + delta*gen["alpha"]
                pg = min(gen["pmax"], pg)
                pg = max(gen["pmin"], pg)
                total_pg += pg
            else
                total_pg += gen["pg_base"]
            end
        end
    end

    return total_pg
end


apply_pg_response!(network, pg_delta::Real) = apply_pg_response!(network, network["response_gens"], pg_delta::Real)

function apply_pg_response!(network, gens::Set{Int}, pg_delta::Real)
    for (i,gen) in network["gen"]
        gen["pg_fixed"] = false
    end

    pg_total = 0.0
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            pg_total += gen["pg"]
        end
    end

    pg_target = pg_total + pg_delta
    #info(LOGGER, "total gen:  $(pg_total)")
    #info(LOGGER, "target gen: $(pg_target)")
    status = 0

    delta_est = 0.0
    while !isapprox(pg_total, pg_target)
        alpha_total = 0.0
        for (i,gen) in network["gen"]
            if gen["gen_status"] != 0 && !gen["pg_fixed"] && gen["index"] in gens
                alpha_total += gen["alpha"]
            end
        end
        #info(LOGGER, "alpha total: $(alpha_total)")

        if isapprox(alpha_total, 0.0) && !isapprox(pg_total, pg_target)
            warn(LOGGER, "insufficient generator response to meet demand, remaining pg $(pg_total - pg_target), remaining alpha $(alpha_total)")
            status = 1
            break
        end

        delta_est += pg_delta/alpha_total
        #info(LOGGER, "detla: $(delta_est)")

        for (i,gen) in network["gen"]
            if gen["gen_status"] != 0 && gen["index"] in gens
                pg_cont = gen["pg_base"] + delta_est*gen["alpha"]

                if pg_cont <= gen["pmin"]
                    gen["pg"] = gen["pmin"]
                    if !gen["pg_fixed"]
                        gen["pg_fixed"] = true
                        #info(LOGGER, "gen $(i) hit lb $(gen["pmin"]) with target value of $(pg_cont)")
                    end
                elseif pg_cont >= gen["pmax"]
                    gen["pg"] = gen["pmax"]
                    if !gen["pg_fixed"]
                        gen["pg_fixed"] = true
                        #info(LOGGER, "gen $(i) hit ub $(gen["pmax"]) with target value of $(pg_cont)")
                    end
                else
                    gen["pg"] = pg_cont
                end
            end
        end

        pg_total = 0.0
        for (i,gen) in network["gen"]
            if gen["gen_status"] != 0
                pg_total += gen["pg"]
            end
        end

        #pg_comp = comp_pg_response_total(network, delta=delta_est)
        #info(LOGGER, "detla: $(delta_est)")
        #info(LOGGER, "total gen comp $(pg_comp) - gen inc $(pg_total)")
        #info(LOGGER, "total gen $(pg_total) - target gen $(pg_target)")

        pg_delta = pg_target - pg_total
    end

    alpha_final = 0.0
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0 && !gen["pg_fixed"] && gen["index"] in gens
            alpha_final += gen["alpha"]
        end
    end
    if isapprox(alpha_final, 0.0)
        warn(LOGGER, "no remaining alpha for generator response (final alpha value $(alpha_final))")
        debug(LOGGER, "delta $(delta_est)")
    end

    network["delta"] = delta_est
    return status
end


"fixes solution degeneracy issues when qg is a free variable, as is the case in PowerFlow"
function correct_qg!(network, solution; bus_gens=gens_by_bus(network))
    for (i,gens) in bus_gens
        if length(gens) > 1
            gen_ids = [gen["index"] for gen in gens]
            qgs = [solution["gen"]["$(j)"]["qg"] for j in gen_ids]
            if !isapprox(abs(sum(qgs)), sum(abs.(qgs)))
                #info(LOGGER, "$(i) - $(gen_ids) - $(qgs) - output requires correction!")
                qg_total = sum(qgs)

                qg_remaining = sum(qgs)
                qg_assignment = Dict(j => 0.0 for j in gen_ids)
                for (i,gen) in enumerate(gens)
                    gen_qg = qg_remaining
                    gen_qg = max(gen_qg, gen["qmin"])
                    gen_qg = min(gen_qg, gen["qmax"])
                    qg_assignment[gen["index"]] = gen_qg
                    qg_remaining = qg_remaining - gen_qg
                    if i == length(gens) && abs(qg_remaining) > 0.0
                        qg_assignment[gen["index"]] = gen_qg + qg_remaining
                    end
                end
                #info(LOGGER, "$(qg_assignment)")
                for (j,qg) in qg_assignment
                    solution["gen"]["$(j)"]["qg"] = qg
                end

                sol_qg_total = sum(solution["gen"]["$(j)"]["qg"] for j in gen_ids)
                #info(LOGGER, "$(qg_total) - $(sol_qg_total)")
                @assert isapprox(qg_total, sol_qg_total)

                #info(LOGGER, "updated to $([solution["gen"]["$(j)"]["qg"] for j in gen_ids])")
            end
        end
    end
end


""
function solution_second_stage!(pm::GenericPowerModel, sol::Dict{String,Any})
    #start_time = time()
    PowerModels.add_setpoint_bus_voltage!(sol, pm)
    #Memento.info(LOGGER, "voltage solution time: $(time() - start_time)")

    #start_time = time()
    PowerModels.add_setpoint_generator_power!(sol, pm)
    #Memento.info(LOGGER, "generator solution time: $(time() - start_time)")

    #start_time = time()
    PowerModels.add_setpoint_branch_flow!(sol, pm)
    #Memento.info(LOGGER, "branch solution time: $(time() - start_time)")

    sol["delta"] = JuMP.value(var(pm, :delta))

    #start_time = time()
    PowerModels.add_setpoint_fixed!(sol, pm, "shunt", "bs", default_value = (item) -> item["bs"])
    #Memento.info(LOGGER, "shunt solution time: $(time() - start_time)")

    #PowerModels.print_summary(sol)
end


# NOTE this does not appear to help, it seems that PF requires bus switching is always true.
function run_fixpoint_pf_soft!(network, pg_lost, model_constructor, solver; iteration_limit=typemax(Int64))
    time_start = time()

    #delta = apply_pg_response!(network, pg_lost)
    #delta = 0.0

    #info(LOGGER, "pg lost: $(pg_lost)")
    #info(LOGGER, "delta: $(network["delta"])")
    #info(LOGGER, "pre-solve time: $(time() - time_start)")

    base_solution = extract_solution(network)
    base_solution["delta"] = network["delta"]
    result = Dict(
        "termination_status" => LOCALLY_SOLVED,
        "solution" => base_solution
    )

    gen_idx_max = maximum(gen["index"] for (i,gen) in network["gen"])
    gen_idx = trunc(Int, 10^ceil(log10(gen_idx_max)))

    bus_gens = gens_by_bus(network)


    iteration = 1
    vm_fixed = true
    while vm_fixed && iteration < iteration_limit
        info(LOGGER, "pf soft fixpoint iteration: $iteration")

        time_start = time()
        result = run_pf_soft_rect(network, model_constructor, solver)
        info(LOGGER, "solve pf time: $(time() - time_start)")

        if result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
            correct_qg!(network, result["solution"], bus_gens=bus_gens)
            PowerModels.update_data!(network, result["solution"])
        else
            warn(LOGGER, "solve issue with run_pf_soft_rect, $(result["termination_status"])")
            break
        end

        vm_fixed = false
        for (i,bus) in network["bus"]
            if bus["vm"] < bus["vmin"] - vm_bound_tol || bus["vm"] > bus["vmax"] + vm_bound_tol
                active_gens = 0
                if length(bus_gens[i]) > 0
                    active_gens = sum(gen["gen_status"] != 0 for gen in bus_gens[i])
                end
                @assert(active_gens == 0)

                bus["vm_fixed"] = true
                current_vm = bus["vm"]
                if bus["vm"] < bus["vmin"]
                    bus["vm_base"] = bus["vmin"]
                    bus["vm"] = bus["vmin"]
                end
                if bus["vm"] > bus["vmax"]
                    bus["vm_base"] = bus["vmax"]
                    bus["vm"] = bus["vmax"]
                end

                warn(LOGGER, "bus $(i) voltage out of bounds $(current_vm) -> $(bus["vm"]), adding virtual generator $(gen_idx)")
                gen_virtual = deepcopy(gen_default)
                gen_virtual["index"] = gen_idx
                gen_virtual["gen_bus"] = bus["index"]
                push!(bus_gens[i], gen_virtual)

                @assert(!haskey(network["gen"], "$(gen_idx)"))
                network["gen"]["$(gen_idx)"] = gen_virtual
                gen_idx += 1

                vm_fixed = true
            end
        end
        iteration += 1
    end

    if iteration >= iteration_limit
        warn(LOGGER, "hit iteration limit")
    end

    for (i,gen) in collect(network["gen"])
        if haskey(gen, "virtual") && gen["virtual"]
            delete!(network["gen"], i)
            delete!(result["solution"]["gen"], i)
        end
    end

    return result
end



""
function run_pf_soft_rect(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_pf_soft_rect; solution_builder = solution_second_stage!, kwargs...)
end

""
function post_pf_soft_rect(pm::GenericPowerModel)
    PowerModels.variable_voltage(pm, bounded=false)
    PowerModels.variable_active_generation(pm, bounded=false)
    PowerModels.variable_reactive_generation(pm, bounded=false)

    delta = ref(pm, :delta)
    var(pm)[:delta] = @variable(pm.model, base_name="delta", start=0.0)
    #@constraint(pm.model, var(pm, :delta) == 0.0)

    var(pm)[:p_slack] = @variable(pm.model, p_slack, base_name="p_slack", start=0.0)

    vr = var(pm, :vr)
    vi = var(pm, :vi)

    PowerModels.constraint_model_voltage(pm)

    for (i,bus) in ref(pm, :bus)
        if bus["vm_fixed"]
            @constraint(pm.model, vr[i]^2 + vi[i]^2 == bus["vm_base"]^2)
        end
    end

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end
    #Memento.info(LOGGER, "misc constraints time: $(time() - start_time)")


    start_time = time()
    p = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    q = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    for (i,branch) in ref(pm, :branch)
        #PowerModels.constraint_ohms_yt_from(pm, i)
        #PowerModels.constraint_ohms_yt_to(pm, i)

        f_bus_id = branch["f_bus"]
        t_bus_id = branch["t_bus"]
        f_idx = (i, f_bus_id, t_bus_id)
        t_idx = (i, t_bus_id, f_bus_id)

        f_bus = ref(pm, :bus, f_bus_id)
        t_bus = ref(pm, :bus, t_bus_id)

        #g, b = PowerModels.calc_branch_y(branch)
        #tr, ti = PowerModels.calc_branch_t(branch)
        g = branch["g"]
        b = branch["b"]
        tr = branch["tr"]
        ti = branch["ti"]

        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        #p_fr  = var(pm, :p, f_idx)
        #q_fr  = var(pm, :q, f_idx)
        #p_to  = var(pm, :p, t_idx)
        #q_to  = var(pm, :q, t_idx)
        vr_fr = vr[f_bus_id] #var(pm, :vr, f_bus_id)
        vr_to = vr[t_bus_id] #var(pm, :vr, t_bus_id)
        vi_fr = vi[f_bus_id] #var(pm, :vi, f_bus_id)
        vi_to = vi[t_bus_id] #var(pm, :vi, t_bus_id)


        if branch["transformer"]
            p[f_idx] = (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        else
            p[f_idx] = (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end
        p[t_idx] = (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))

    end
    #Memento.info(LOGGER, "flow expr time: $(time() - start_time)")


    start_time = time()
    pg = Dict{Int,Any}()
    qg = Dict{Int,Any}()
    for (i,gen) in ref(pm, :gen)
        #=
        if gen["pg_fixed"]
            pg[i] = gen["pg"]
            @constraint(pm.model, var(pm, :pg, i) == gen["pg"])
        else
            if i in ref(pm, :response_gens)
                pg[i] = gen["pg_base"] + gen["alpha"]*delta
                @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"] + gen["alpha"]*delta)
            else
                pg[i] = gen["pg_base"]
                @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"])
            end
        end
        =#
        pg[i] = gen["pg_base"]
        @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"])

        #=
        if gen["qg_fixed"]
            qg[i] = gen["qg"]
            @constraint(pm.model, var(pm, :qg, i) == gen["qg"])
        else
            qg[i] = var(pm, :qg, i)
        end
        =#
        qg[i] = var(pm, :qg, i)
    end
    #Memento.info(LOGGER, "gen expr time: $(time() - start_time)")


    start_time = time()
    for (i,bus) in ref(pm, :bus)
        #PowerModels.constraint_power_balance_shunt(pm, i)

        bus_arcs = ref(pm, :bus_arcs, i)
        bus_arcs_dc = ref(pm, :bus_arcs_dc, i)
        bus_gens = ref(pm, :bus_gens, i)
        bus_loads = ref(pm, :bus_loads, i)
        bus_shunts = ref(pm, :bus_shunts, i)

        bus_pd = Dict(k => ref(pm, :load, k, "pd") for k in bus_loads)
        bus_qd = Dict(k => ref(pm, :load, k, "qd") for k in bus_loads)

        bus_gs = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts)
        bus_bs = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts)

        #p = var(pm, :p)
        #q = var(pm, :q)
        #pg = var(pm, :pg)
        #qg = var(pm, :qg)

        @constraint(pm.model, p_slack + sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr[i]^2 + vi[i]^2))
        @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr[i]^2 + vi[i]^2))
    end
    #Memento.info(LOGGER, "power balance constraint time: $(time() - start_time)")
end








function run_fixpoint_pf_v2_3!(network, pg_lost, model_constructor, solver; iteration_limit=typemax(Int64))
    time_start = time()

    network_backup = deepcopy(network)

    response_status = apply_pg_response!(network, pg_lost)

    debug(LOGGER, "pg lost: $(pg_lost)")
    debug(LOGGER, "delta: $(network["delta"])")
    debug(LOGGER, "status: $(response_status)")
    debug(LOGGER, "pre-solve time: $(time() - time_start)")


    base_solution = extract_solution(network)
    base_solution["delta"] = network["delta"]
    final_result = Dict(
        "termination_status" => LOCALLY_SOLVED,
        "solution" => base_solution
    )

    bus_gens = gens_by_bus(network)

    active_response_gens = []
    for i in network["response_gens"]
        gen = network["gen"]["$(i)"]
        if gen["gen_status"] != 0
            push!(active_response_gens, gen)
        end
    end

    #network["delta_start"] = network["delta"]
    for (i,gen) in network["gen"]
        gen["qg_fixed"] = false
        gen["pg_start"] = gen["pg"]
        if isapprox(gen["qmin"],gen["qmax"])
            gen["qg_fixed"] = true
            gen["qg"] = gen["qmin"]
        end
        gen["qg_start"] = gen["qg"]
    end

    for (i,bus) in network["bus"]
        active_gens = [gen for gen in bus_gens[i] if !gen["qg_fixed"]]
        if length(active_gens) == 0
            bus["vm_fixed"] = false
        else
            bus["vm_fixed"] = true
        end
        #bus["vm_start"] = bus["vm"]
        #bus["va_start"] = bus["va"]
        bus["vr_start"] = bus["vm"]*cos(bus["va"])
        bus["vi_start"] = bus["vm"]*sin(bus["va"])
    end

    pf_fixed_all = all(gen["pg_fixed"] for gen in active_response_gens)

    cont_pf_failed = false
    time_start = time()
    #result = run_fixed_pf_nbf_rect(network, model_constructor, solver)
    if !pf_fixed_all
        result = run_fixed_pf_nbf_rect2(network, model_constructor, solver)
    else
        result = run_fixed_pf_nbf_rect2_ds(network, model_constructor, solver)
    end
    debug(LOGGER, "pf solve time: $(time() - time_start)")
    if result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
        correct_qg!(network, result["solution"], bus_gens=bus_gens)
        PowerModels.update_data!(network, result["solution"])
        final_result = result
    else
        warn(LOGGER, "$(network["cont_label"]) contingency pf solver FAILED with status $(result["termination_status"]) on iteration 0")
        cont_pf_failed = true
    end



    pg_switched = true
    qg_switched = true
    vm_switched = true

    iteration = 1
    deltas = [result["solution"]["delta"]]
    while (pg_switched || qg_switched || vm_switched) && !cont_pf_failed && iteration <= iteration_limit
        debug(LOGGER, "obj: $(result["objective"])")
        debug(LOGGER, "delta: $(result["solution"]["delta"])")
        pg_switched = false
        qg_switched = false
        vm_switched = false

        for (i,gen) in network["gen"]
            if gen["index"] in network["response_gens"]
                pg = gen["pg_base"] + network["delta"]*gen["alpha"]
            else
                pg = gen["pg_base"]
            end

            if gen["pg_fixed"]
                if !isapprox(gen["pmax"], gen["pmin"]) && pg < gen["pmax"] && pg > gen["pmin"]
                    gen["pg"] = pg
                    gen["pg_fixed"] = false
                    pg_switched = true
                    #info(LOGGER, "unfix pg on gen $(i)")
                end
            else
                if pg >= gen["pmax"]
                    gen["pg"] = gen["pmax"]
                    gen["pg_fixed"] = true
                    pg_switched = true
                    #info(LOGGER, "fix pg to ub on gen $(i)")
                elseif gen["pg"] <= gen["pmin"]
                    gen["pg"] = gen["pmin"]
                    gen["pg_fixed"] = true
                    pg_switched = true
                    #info(LOGGER, "fix pg to lb on gen $(i)")
                end
            end
        end

        for (i,bus) in network["bus"]
            if length(bus_gens[i]) > 0
                qg = sum(gen["qg"] for gen in bus_gens[i])
                qmin = sum(gen["qmin"] for gen in bus_gens[i])
                qmax = sum(gen["qmax"] for gen in bus_gens[i])

                if isapprox(qmin,qmax)
                    @assert !bus["vm_fixed"]
                    for gen in bus_gens[i]
                        @assert gen["qg_fixed"]
                        @assert isapprox(gen["qg"],gen["qmin"])
                    end
                elseif bus["vm_fixed"]
                    if qg >= qmax
                        bus["vm_fixed"] = false
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg"] = gen["qmax"]
                            gen["qg_fixed"] = true
                        end
                        #info(LOGGER, "fix qg to ub on bus $(i)")
                    end

                    if qg <= qmin
                        bus["vm_fixed"] = false
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg"] = gen["qmin"]
                            gen["qg_fixed"] = true
                        end
                        #info(LOGGER, "fix qg to lb on bus $(i)")
                    end
                else
                    if qg < qmax && qg > qmin
                        bus["vm_fixed"] = true

                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg_fixed"] = false
                            gen["qg_start"] = gen["qg"]
                        end
                        #info(LOGGER, "fix vm to on bus $(i)")
                    end
                    if qg >= qmax && bus["vm"] > bus["vm_base"]
                        bus["vm_fixed"] = true
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg_fixed"] = false
                        end
                        #info(LOGGER, "fix vm to on bus $(i)")
                    end
                    if qg <= qmin && bus["vm"] < bus["vm_base"]
                        bus["vm_fixed"] = true
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg_fixed"] = false
                        end
                        #info(LOGGER, "fix vm to on bus $(i)")
                    end
                end

                #=
                leads to solver infeasiblity
                for gen in bus_gens[i]
                    if isapprox(gen["qmin"],gen["qmax"])
                        gen["qg_fixed"] = true
                        gen["qg"] = gen["qmin"]
                        gen["qg_start"] = gen["qg"]
                    end
                end
                =#
            end
        end


        for (i,gen) in network["gen"]
            gen["pg_start"] = gen["pg"]
            gen["qg_start"] = gen["qg"]
        end

        for (i,bus) in network["bus"]
            bus["vm_start"] = bus["vm"]
            bus["va_start"] = bus["va"]
        end


        if pg_switched || qg_switched || vm_switched
            debug(LOGGER, "bus or gen switched: $iteration")
            time_start = time()

            pf_fixed_all = all(gen["pg_fixed"] for gen in active_response_gens)

            #result = run_fixed_pf_nbf_rect(network, model_constructor, solver)
            #result = run_fixed_pf_nbf_rect2(network, model_constructor, solver)
            if !pf_fixed_all
                result = run_fixed_pf_nbf_rect2(network, model_constructor, solver)
            else
                result = run_fixed_pf_nbf_rect2_ds(network, model_constructor, solver)
            end
            debug(LOGGER, "pf solve time: $(time() - time_start)")
            if result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
                correct_qg!(network, result["solution"], bus_gens=bus_gens)
                PowerModels.update_data!(network, result["solution"])
                final_result = result
            else
                warn(LOGGER, "$(network["cont_label"]) contingency pf solver FAILED with status $(result["termination_status"]) on iteration 0")
                break
            end

            push!(deltas, result["solution"]["delta"])
            iteration += 1
            if iteration >= iteration_limit
                warn(LOGGER, "hit iteration limit")
            end
            if length(deltas) > 3 && isapprox(deltas[end-2], deltas[end])
                warn(LOGGER, "cycle detected, stopping")
                break
            end
        end
    end

    vm_bound_vio = false
    for (i,bus) in network["bus"]
        bus_sol = final_result["solution"]["bus"][i]
        if bus_sol["vm"] - vm_bound_tol >= bus["vmax"] || bus_sol["vm"] + vm_bound_tol <= bus["vmin"]
            vm_bound_vio = true
            warn(LOGGER, "$(network["cont_label"]) vm bound out of range on bus $(i): $(bus["vmin"]) - $(bus_sol["vm"]) - $(bus["vmax"])")
        end
    end

    qg_bound_vio = false
    for (i,gen) in network["gen"]
        gen_sol = final_result["solution"]["gen"][i]
        if gen_sol["qg"] - qg_bound_tol >= gen["qmax"] || gen_sol["qg"] + qg_bound_tol <= gen["qmin"]
            qg_bound_vio = true
            warn(LOGGER, "$(network["cont_label"]) qg bound out of range on gen $(i): $(gen["qmin"]) - $(gen_sol["qg"]) - $(gen["qmax"])")
        end
    end

    if vm_bound_vio || qg_bound_vio
        warn(LOGGER, "$(network["cont_label"]) running voltage profile correction")
        result = run_fixpoint_pf_v5!(network_backup, pg_lost, model_constructor, solver, iteration_limit=iteration_limit)
        if result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
            PowerModels.update_data!(network, result["solution"])
            final_result = result
        else
            warn(LOGGER, "$(network["cont_label"]) voltage profile correction solver FAILED with status $(result["termination_status"])")
        end
    end

    return final_result
end




""
function run_fixed_pf_nbf_rect2(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_fixed_pf_nbf_rect2; solution_builder = solution_second_stage!, kwargs...)
end

""
function post_fixed_pf_nbf_rect2(pm::GenericPowerModel)
    start_time = time()
    PowerModels.variable_voltage(pm, bounded=false)
    PowerModels.variable_active_generation(pm, bounded=false)
    PowerModels.variable_reactive_generation(pm, bounded=false)
    #PowerModels.variable_branch_flow(pm, bounded=false)
    #PowerModels.variable_dcline_flow(pm, bounded=false)

    #PowerModels.variable_branch_flow(pm, bounded=false)

    # TODO set bounds bounds on alpha and total gen capacity
    var(pm)[:delta] = @variable(pm.model, delta, base_name="delta", start=0.0)
    #var(pm)[:delta] = @variable(pm.model, delta, base_name="delta", start=ref(pm, :delta_start))
    #Memento.info(LOGGER, "post variable time: $(time() - start_time)")

    start_time = time()

    vr = var(pm, :vr)
    vi = var(pm, :vi)

    PowerModels.constraint_model_voltage(pm)

    for (i,bus) in ref(pm, :bus)
        if bus["vm_fixed"]
            @constraint(pm.model, vr[i]^2 + vi[i]^2 == bus["vm_base"]^2)
        end
    end

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end
    #Memento.info(LOGGER, "misc constraints time: $(time() - start_time)")


    start_time = time()
    p = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    q = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    for (i,branch) in ref(pm, :branch)
        #PowerModels.constraint_ohms_yt_from(pm, i)
        #PowerModels.constraint_ohms_yt_to(pm, i)

        f_bus_id = branch["f_bus"]
        t_bus_id = branch["t_bus"]
        f_idx = (i, f_bus_id, t_bus_id)
        t_idx = (i, t_bus_id, f_bus_id)

        f_bus = ref(pm, :bus, f_bus_id)
        t_bus = ref(pm, :bus, t_bus_id)

        #g, b = PowerModels.calc_branch_y(branch)
        #tr, ti = PowerModels.calc_branch_t(branch)
        g = branch["g"]
        b = branch["b"]
        tr = branch["tr"]
        ti = branch["ti"]

        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        #p_fr  = var(pm, :p, f_idx)
        #q_fr  = var(pm, :q, f_idx)
        #p_to  = var(pm, :p, t_idx)
        #q_to  = var(pm, :q, t_idx)
        vr_fr = vr[f_bus_id] #var(pm, :vr, f_bus_id)
        vr_to = vr[t_bus_id] #var(pm, :vr, t_bus_id)
        vi_fr = vi[f_bus_id] #var(pm, :vi, f_bus_id)
        vi_to = vi[t_bus_id] #var(pm, :vi, t_bus_id)


        if branch["transformer"]
            p[f_idx] = (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        else
            p[f_idx] = (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end
        p[t_idx] = (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))

    end
    #Memento.info(LOGGER, "flow expr time: $(time() - start_time)")


    start_time = time()
    pg = Dict{Int,Any}()
    qg = Dict{Int,Any}()
    for (i,gen) in ref(pm, :gen)
        if gen["pg_fixed"]
            pg[i] = gen["pg"]
            @constraint(pm.model, var(pm, :pg, i) == gen["pg"])
        else
            if i in ref(pm, :response_gens)
                pg[i] = gen["pg_base"] + gen["alpha"]*delta
                @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"] + gen["alpha"]*delta)
            else
                pg[i] = gen["pg_base"]
                @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"])
            end
        end

        if gen["qg_fixed"]
            qg[i] = gen["qg"]
            @constraint(pm.model, var(pm, :qg, i) == gen["qg"])
        else
            qg[i] = var(pm, :qg, i)
        end
    end
    #Memento.info(LOGGER, "gen expr time: $(time() - start_time)")


    start_time = time()
    for (i,bus) in ref(pm, :bus)
        #PowerModels.constraint_power_balance_shunt(pm, i)

        bus_arcs = ref(pm, :bus_arcs, i)
        bus_arcs_dc = ref(pm, :bus_arcs_dc, i)
        bus_gens = ref(pm, :bus_gens, i)
        bus_loads = ref(pm, :bus_loads, i)
        bus_shunts = ref(pm, :bus_shunts, i)

        bus_pd = Dict(k => ref(pm, :load, k, "pd") for k in bus_loads)
        bus_qd = Dict(k => ref(pm, :load, k, "qd") for k in bus_loads)

        bus_gs = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts)
        bus_bs = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts)

        #p = var(pm, :p)
        #q = var(pm, :q)
        #pg = var(pm, :pg)
        #qg = var(pm, :qg)

        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr[i]^2 + vi[i]^2))
        @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr[i]^2 + vi[i]^2))
    end
    #Memento.info(LOGGER, "power balance constraint time: $(time() - start_time)")
end


"a variant of fixed_pf_nbf_rect2 with a distributed active power slack"
function run_fixed_pf_nbf_rect2_ds(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_fixed_pf_nbf_rect2_ds; solution_builder = solution_second_stage!, kwargs...)
end

""
function post_fixed_pf_nbf_rect2_ds(pm::GenericPowerModel)
    PowerModels.variable_voltage(pm, bounded=false)
    PowerModels.variable_active_generation(pm, bounded=false)
    PowerModels.variable_reactive_generation(pm, bounded=false)
    #PowerModels.variable_branch_flow(pm, bounded=false)
    #PowerModels.variable_dcline_flow(pm, bounded=false)

    #PowerModels.variable_branch_flow(pm, bounded=false)

    delta = ref(pm, :delta)
    var(pm)[:delta] = @variable(pm.model, base_name="delta", start=0.0)
    @constraint(pm.model, var(pm, :delta) == delta)

    active_response_gens = intersect(ids(pm, :gen), ref(pm, :response_gens))
    # for i in active_response_gens
    #     gen = ref(pm, :gen, i)
    #     println("$(i) $(gen["pg_fixed"]) $(gen["pmin"]) $(gen["pg_base"]) $(gen["pg"]) $(gen["pmax"])")
    # end

    if !all(ref(pm, :gen, i)["pg_fixed"] for i in active_response_gens)
        Memento.error(LOGGER, "fixed_pf_nbf_rect2_ds model requires all response_gens have pg_fixed set to true")
    end
    var(pm)[:p_slack] = @variable(pm.model, p_slack, base_name="p_slack", start=0.0)

    vr = var(pm, :vr)
    vi = var(pm, :vi)

    PowerModels.constraint_model_voltage(pm)

    for (i,bus) in ref(pm, :bus)
        if bus["vm_fixed"]
            @constraint(pm.model, vr[i]^2 + vi[i]^2 == bus["vm_base"]^2)
        end
    end

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end
    #Memento.info(LOGGER, "misc constraints time: $(time() - start_time)")


    start_time = time()
    p = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    q = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    for (i,branch) in ref(pm, :branch)
        #PowerModels.constraint_ohms_yt_from(pm, i)
        #PowerModels.constraint_ohms_yt_to(pm, i)

        f_bus_id = branch["f_bus"]
        t_bus_id = branch["t_bus"]
        f_idx = (i, f_bus_id, t_bus_id)
        t_idx = (i, t_bus_id, f_bus_id)

        f_bus = ref(pm, :bus, f_bus_id)
        t_bus = ref(pm, :bus, t_bus_id)

        #g, b = PowerModels.calc_branch_y(branch)
        #tr, ti = PowerModels.calc_branch_t(branch)
        g = branch["g"]
        b = branch["b"]
        tr = branch["tr"]
        ti = branch["ti"]

        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        #p_fr  = var(pm, :p, f_idx)
        #q_fr  = var(pm, :q, f_idx)
        #p_to  = var(pm, :p, t_idx)
        #q_to  = var(pm, :q, t_idx)
        vr_fr = vr[f_bus_id] #var(pm, :vr, f_bus_id)
        vr_to = vr[t_bus_id] #var(pm, :vr, t_bus_id)
        vi_fr = vi[f_bus_id] #var(pm, :vi, f_bus_id)
        vi_to = vi[t_bus_id] #var(pm, :vi, t_bus_id)


        if branch["transformer"]
            p[f_idx] = (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        else
            p[f_idx] = (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end
        p[t_idx] = (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))

    end
    #Memento.info(LOGGER, "flow expr time: $(time() - start_time)")


    start_time = time()
    pg = Dict{Int,Any}()
    qg = Dict{Int,Any}()
    for (i,gen) in ref(pm, :gen)
        if gen["pg_fixed"]
            pg[i] = gen["pg"]
            @constraint(pm.model, var(pm, :pg, i) == gen["pg"])
        else
            if i in ref(pm, :response_gens)
                pg[i] = gen["pg_base"] + gen["alpha"]*delta
                @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"] + gen["alpha"]*delta)
            else
                pg[i] = gen["pg_base"]
                @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"])
            end
        end

        if gen["qg_fixed"]
            qg[i] = gen["qg"]
            @constraint(pm.model, var(pm, :qg, i) == gen["qg"])
        else
            qg[i] = var(pm, :qg, i)
        end
    end
    #Memento.info(LOGGER, "gen expr time: $(time() - start_time)")


    start_time = time()
    for (i,bus) in ref(pm, :bus)
        #PowerModels.constraint_power_balance_shunt(pm, i)

        bus_arcs = ref(pm, :bus_arcs, i)
        bus_arcs_dc = ref(pm, :bus_arcs_dc, i)
        bus_gens = ref(pm, :bus_gens, i)
        bus_loads = ref(pm, :bus_loads, i)
        bus_shunts = ref(pm, :bus_shunts, i)

        bus_pd = Dict(k => ref(pm, :load, k, "pd") for k in bus_loads)
        bus_qd = Dict(k => ref(pm, :load, k, "qd") for k in bus_loads)

        bus_gs = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts)
        bus_bs = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts)

        #p = var(pm, :p)
        #q = var(pm, :q)
        #pg = var(pm, :pg)
        #qg = var(pm, :qg)

        @constraint(pm.model, p_slack + sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr[i]^2 + vi[i]^2))
        @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr[i]^2 + vi[i]^2))
    end
    #Memento.info(LOGGER, "power balance constraint time: $(time() - start_time)")
end

