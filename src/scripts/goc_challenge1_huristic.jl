#!/usr/bin/env julia

#start_pkg = time()

include("distributed.jl")

#include("second-stage-fp.jl")
include("second-stage-soft-fp.jl")

using Memento
const LOGGER = Memento.getlogger(PowerModelsSecurityConstrained)

using JuMP
using PowerModels

using Ipopt
#using Gurobi
using Cbc

#@everywhere Memento.config("debug")
#@everywhere setlevel!(LOGGER, "debug")

#println("package load time: $(time() - start_pkg)")

#println("script startup time: $(time() - start_init)")


function compute_solution1(con_file::String, inl_file::String, raw_file::String, rop_file::String, time_limit::Int, scoring_method::Int, network_model::String; output_dir::String="", scenario_id::String="none")
    time_start = time()
    info(LOGGER, "time remaining: $(time_limit)")

    goc_data = parse_goc_files(con_file, inl_file, raw_file, rop_file, scenario_id=scenario_id)
    network = build_pm_model(goc_data)
    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    correct_network_solution!(network)
    write_solution1(network, output_dir=output_dir)

    result = Dict(
        "termination_status" => LOCALLY_SOLVED,
        "solution" => extract_solution(network)
    )

    time_data = time() - time_start


    ###### Preparations for Solver ######

    time_solve_start = time()
    nlp_solver = with_optimizer(Ipopt.Optimizer, tol=1e-6, mu_init=1e1)
    nlp_solver_relaxed = with_optimizer(Ipopt.Optimizer, tol=1e-6, mu_init=1e1)
    #qp_solver = with_optimizer(Gurobi.Optimizer, OptimalityTol=1e-6, Method=2, Crossover=0)
    #qp_solver_relaxed = with_optimizer(Gurobi.Optimizer, OptimalityTol=1e-6, Method=2, Crossover=0, BarConvTol=5e-3)
    qp_solver = with_optimizer(Cbc.Optimizer)
    qp_solver_relaxed = qp_solver
    lp_solver = qp_solver

    time_filter = 0.0

    ###### DC OPF Solve ######

    time_dc_opf_start = time()

    network_apo = deepcopy(network)

    result = run_opf_cheap_dc(network_apo, DCPPowerModel, lp_solver)
    if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
        warn(LOGGER, "base case DC-OPF solve failed with status $(result["termination_status"]), try with relaxed convergence tolerance")

        result = run_opf_cheap_dc(network_apo, DCPPowerModel, qp_solver_relaxed)
        if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
            warn(LOGGER, "relaxed base case DC-OPF solve failed with status $(result["termination_status"])")
            result = Dict(
                "termination_status" => LOCALLY_SOLVED,
                "solution" => extract_solution(network)
            )
        else
            warn(LOGGER, "relaxed base case DC-OPF solve status $(result["termination_status"])")
        end
    end

    update_active_power_data!(network_apo, result["solution"])
    update_active_power_data!(network, result["solution"])


    write_solution1(network_apo, output_dir=output_dir, solution_file="solution1_apo.txt")

    correct_network_solution!(network)
    write_solution1(network, output_dir=output_dir)

    for (i,bus) in network["bus"]
        bus["vm"] = (bus["vmax"] + bus["vmin"])/2.0 + 0.04
        bus["vm_start"] = bus["vm"]
        bus["va_start"] = bus["va"]

        bus["vr_start"] = bus["vm"]*cos(bus["va"])
        bus["vi_start"] = bus["vm"]*sin(bus["va"])
    end
    for (i,gen) in network["gen"]
        gen["qg_start"] = 0.0
        gen["pg_start"] = gen["pg"]
    end

    gen_cost = PowerModels.calc_gen_cost(network)
    info(LOGGER, "generation cost: $(gen_cost)")

    time_dc_opf = time() - time_dc_opf_start
    info(LOGGER, "dc opf phase time: $(time_dc_opf)")



    ###### AC OPF Solve ######

    time_ac_opf_start = time()

    #result = run_opf_cheap(network, ACPPowerModel, nlp_solver)
    # if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED || result["termination_status"] == :Suboptimal)
    #     error(LOGGER, "voltage profile solver failed")
    # end

    deactivate_rate_a!(network)
    activate_rate_a_violations!(network)

    line_flow_vio = true
    while line_flow_vio

        result = run_opf_pg_pf_rect_5(network, ACRPowerModel, nlp_solver)
        if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED)
            warn(LOGGER, "base case AC-OPF solve failed with status $(result["termination_status"]), try with relaxed convergence tolerance")
            break
            # result = run_opf_pg_pf_rect_5(network, ACRPowerModel, nlp_solver_relaxed)
            # if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
            #     warn(LOGGER, "relaxed base case AC-OPF solve failed with status $(result["termination_status"])")
            #     break
            # else
            #     warn(LOGGER, "relaxed base case AC-OPF solve status $(result["termination_status"])")
            # end
        end
        PowerModels.update_data!(network, result["solution"])

        for (i,bus) in network["bus"]
            bus["vm_start"] = bus["vm"]
            bus["va_start"] = bus["va"]

            bus["vr_start"] = bus["vm"]*cos(bus["va"])
            bus["vi_start"] = bus["vm"]*sin(bus["va"])
        end
        for (i,gen) in network["gen"]
            gen["qg_start"] = gen["qg"]
            gen["pg_start"] = gen["pg"]
        end
        line_flow_vio = activate_rate_a_violations!(network)
    end
    #activate_rate_a!(network)

    balance = compute_power_balance_deltas!(network)
    info(LOGGER, "power balance cost: $(balance)")


    #PowerModels.update_data!(network, result["solution"])
    correct_network_solution!(network)
    write_solution1(network, output_dir=output_dir)

    time_ac_opf = time() - time_ac_opf_start
    info(LOGGER, "ac opf phase time: $(time_ac_opf)")


    network["delta"] = 0.0
    for (i,bus) in network["bus"]
        bus["va_start"] = bus["va"]
        bus["vm_start"] = bus["vm"]
    end

    for (i,gen) in network["gen"]
        gen["pg_base"] = gen["pg"]
        gen["pg_start"] = gen["pg"]
        gen["qg_start"] = gen["qg"]
    end

    time_solve = time() - time_solve_start
    info(LOGGER, "total solve time: $(time_solve)")


    ###### Memory Check######
    scopf_comp = true
    inverse_size = compute_inverse_size(network)
    workers = Distributed.workers()

    if 2*length(workers)*inverse_size > Sys.total_memory()
        scopf_comp = false
        gb_scale = 2^30
        inverse_mem = trunc(inverse_size/gb_scale, digits=2)
        required_mem = trunc((2*length(workers)*inverse_size)/gb_scale, digits=2)
        system_mem = trunc(Sys.total_memory()/gb_scale, digits=2)
        info(LOGGER, "number of workers $(length(workers)), estimated memory per worker $(2*inverse_mem) Gb")
        info(LOGGER, "insufficient memory to run SCOPF cut generation, required $(required_mem) Gb, system $(system_mem) Gb")
    end


    if scopf_comp
        #return
        ###### Setup Workers ######

        time_worker_start = time()
        workers = Distributed.workers()

        info(LOGGER, "start warmup on $(length(workers)) workers")
        worker_futures = []
        for wid in workers
            future = remotecall(load_network_global, wid, con_file, inl_file, raw_file, rop_file, scenario_id)
            push!(worker_futures, future)
        end

        # setup for contigency solve
        gen_cont_total = length(network["gen_contingencies"])
        branch_cont_total = length(network["branch_contingencies"])
        cont_total = gen_cont_total + branch_cont_total
        cont_per_proc = cont_total/length(workers)

        cont_order = contingency_order(network)
        cont_range = []
        for p in 1:length(workers)
            cont_start = trunc(Int, ceil(1+(p-1)*cont_per_proc))
            cont_end = min(cont_total, trunc(Int,ceil(p*cont_per_proc)))
            push!(cont_range, cont_start:cont_end,)
        end

        for (i,rng) in enumerate(cont_range)
            info(LOGGER, "task $(i): $(length(rng)) / $(rng)")
        end
        #pmap(filter_network_global_contingencies, cont_range)
        output_dirs = [output_dir for i in 1:length(workers)]

        info(LOGGER, "waiting for worker warmup to complete: $(time() - time_start)")
        for future in worker_futures
            wait(future)
        end

        time_worker = time() - time_worker_start
        info(LOGGER, "total worker warmup time: $(time_worker)")


        #return
        ###### DC SCOPF Solve ######

        time_remaining = time_limit - (time() - time_start)
        info(LOGGER, "time remaining before cut enumeration: $(time_remaining)")

        time_filter = 0.0
        objective = result["objective"]
        sol = result["solution"]

        iteration = 0
        network_apo["gen_flow_cuts"] = []
        network_apo["branch_flow_cuts"] = []
        solution_file_apo = ["solution1_apo.txt" for p in 1:length(workers)]
        while true
            time_filter_start = time()
            time_start_iteration = time()
            iteration += 1
            info(LOGGER, "cut enumeration iteration: $(iteration)")

            write_active_flow_cuts(network_apo, output_dir=output_dir)
            #cuts = pmap(check_contingencies_branch_flow_remote, cont_range, output_dirs, [iteration for p in 1:length(workers)], [true for p in 1:length(workers)], solution_file_apo)
            #cuts = pmap(check_contingencies_branch_flow_remote_nd, cont_range, output_dirs, [iteration for p in 1:length(workers)], solution_file_apo)
            #cuts = pmap(check_contingencies_branch_flow_remote_nd_first, cont_range, output_dirs, [iteration for p in 1:length(workers)], solution_file_apo)
            cuts = pmap(check_contingencies_branch_flow_remote_nd_first_lazy, cont_range, output_dirs, [iteration for p in 1:length(workers)], solution_file_apo)
            time_filter += time() - time_filter_start

            cuts_found = sum(length(c.gen_cuts)+length(c.branch_cuts) for c in cuts)
            if cuts_found <= 0
                info(LOGGER, "no violated cuts found scopf fixed-point reached")
                network["gen_flow_cuts"] = network_apo["gen_flow_cuts"]
                network["branch_flow_cuts"] = network_apo["branch_flow_cuts"]
                break
            else
                for cut in cuts
                    info(LOGGER, "found $(length(cut.gen_cuts) + length(cut.branch_cuts)) branch flow violations")
                end
            end

            for cut in cuts
                append!(network_apo["gen_flow_cuts"], cut.gen_cuts)
                append!(network_apo["branch_flow_cuts"], cut.branch_cuts)
            end

            info(LOGGER, "active cuts: gen $(length(network_apo["gen_flow_cuts"])), branch $(length(network_apo["branch_flow_cuts"]))")
            #for c in network_apo["gen_contingencies"]
            #    info(LOGGER, "  $(c.label)")
            #end
            #for c in network_apo["branch_contingencies"]
            #    info(LOGGER, "  $(c.label)")
            #end

            time_solve_start = time()
            #result = run_scopf_cuts_dc_soft(network_apo, DCPPowerModel, qp_solver)

            result = run_scopf_cuts_dc_soft_2(network_apo, DCPPowerModel, qp_solver)
            if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
                warn(LOGGER, "scopf solve failed with status $(result["termination_status"])")

                result = run_scopf_cuts_dc_soft_2(network_apo, DCPPowerModel, qp_solver_relaxed)
                if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
                    warn(LOGGER, "relaxed scopf solve failed with status $(result["termination_status"])")
                    break
                else
                    warn(LOGGER, "relaxed scopf solve status $(result["termination_status"])")
                end
            end
            info(LOGGER, "objective: $(result["objective"])")
            time_solve += time() - time_solve_start
            update_active_power_data!(network_apo, result["solution"])
            update_active_power_data!(network, result["solution"])
            write_solution1(network_apo, output_dir=output_dir, solution_file="solution1_apo.txt")

            balance = compute_power_balance_deltas!(network)
            info(LOGGER, "power balance cost: $(balance)")

            gen_cost = PowerModels.calc_gen_cost(network)
            info(LOGGER, "generation cost: $(gen_cost)")

            time_iteration = time() - time_start_iteration
            info(LOGGER, "iteration time: $(time_iteration)")

            time_remaining = time_limit - (time() - time_start)
            info(LOGGER, "time remaining: $(time_remaining)")

            if time_remaining < time_iteration + time_ac_opf
                warn(LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), iteration time $(time_iteration), ac solve time $(time_ac_opf)")
                break
            end
        end


        ###### AC OPF Recovery Solve ######

        # result = run_opf_cheap_target(network, ACPPowerModel, nlp_solver)
        # if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED || result["termination_status"] == :Suboptimal)
        #     error(LOGGER, "voltage profile solver failed")
        # end
        # sol = result["solution"]
        # PowerModels.update_data!(network, sol)

        line_flow_vio = true
        while line_flow_vio
            result = run_opf_pg_pf_rect_5(network, ACRPowerModel, nlp_solver)

            if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED)
                warn(LOGGER, "base case AC polish solve failed with status $(result["termination_status"])")
                break
                # result = run_opf_pg_pf_rect_5(network, ACRPowerModel, nlp_solver_relaxed)
                # if !(result["termination_status"] == OPTIMAL || result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
                #     warn(LOGGER, "relaxed base case AC-OPF solve failed with status $(result["termination_status"])")
                #     break
                # else
                #     warn(LOGGER, "relaxed base case AC-OPF solve status $(result["termination_status"])")
                # end
            end
            PowerModels.update_data!(network, result["solution"])

            for (i,bus) in network["bus"]
                bus["vm_start"] = bus["vm"]
                bus["va_start"] = bus["va"]

                bus["vr_start"] = bus["vm"]*cos(bus["va"])
                bus["vi_start"] = bus["vm"]*sin(bus["va"])
            end
            for (i,gen) in network["gen"]
                gen["qg_start"] = gen["qg"]
                gen["pg_start"] = gen["pg"]
            end
            line_flow_vio = activate_rate_a_violations!(network)
        end

        balance = compute_power_balance_deltas!(network)
        info(LOGGER, "power balance cost: $(balance)")

        gen_cost = PowerModels.calc_gen_cost(network)
        info(LOGGER, "generation cost: $(gen_cost)")

        correct_network_solution!(network)
        write_solution1(network, output_dir=output_dir)
    end

    ###### Results Summary ######

    active_contingencies = Set(cut.cont_label for cut in [network["gen_flow_cuts"]; network["branch_flow_cuts"]])
    network["gen_contingencies_active"] = [cont for cont in network["gen_contingencies"] if cont.label in active_contingencies]
    network["branch_contingencies_active"] = [cont for cont in network["branch_contingencies"] if cont.label in active_contingencies]

    total_cuts = length(network["gen_flow_cuts"]) + length(network["branch_flow_cuts"])

    write_contingencies(network, output_dir=output_dir)

    write_scopf_summary(goc_data.scenario, network, gen_cost, branch_flow_cuts=total_cuts, load_time=time_data, solve_time=time_solve, filter_time=time_filter, total_time = time() - time_start)
end
