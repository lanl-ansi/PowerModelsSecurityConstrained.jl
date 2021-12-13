#!/usr/bin/env julia

include("distributed.jl")
add_procs()

using PowerModels
using Ipopt

using PowerModelsSecurityConstrained

using Memento

using JuMP
using JSON

const _IM = PowerModels._IM
const _PM = PowerModels
const _LOGGER = PowerModelsSecurityConstrained._LOGGER



function code1(con_file::String, json_file::String, raw_file::String, reserved::String, time_limit::Int64, scoring_method::Int64, case_id::String, scenario_id::String; output_dir=".")
    ###### Load Data ######
    time_load_start = time()

    goc_data = parse_c2_files(raw_file, con_file, json_file, case_id=case_id, scenario_id=scenario_id)

    pm_data = build_c2_pm_model(goc_data)
    #print_summary(pm_data)

    correct_c2_solution!(pm_data)
    write_c2_solution(pm_data, output_dir=output_dir)

    time_load = time() - time_load_start


    time_basecase_start = time()

    objective_best = typemin(Float64)

    println("*** Start OPF Solve:")

    normalize_tm_values!(pm_data)
    set_va_start_values!(pm_data)

    nlp_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6)

    result = run_c2_opf_soft(pm_data, ACPPowerModel, nlp_solver)

    if (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == OPTIMAL)
        println("OPF objective: $(result["objective"])")
        update_data!(pm_data, result["solution"])

        correct_c2_solution!(pm_data)
        write_c2_solution(pm_data, output_dir=output_dir)
        objective_best = result["objective"]
    else
        warn(_LOGGER, "run_c2_opf_soft failed with status $(result["termination_status"])")
    end


    ###### Prepare UC Solution ######
    println("")
    println("*** Start UC-OPF Solve:")

    #println([gen["gen_status"] for (i,gen) in pm_data["gen"]])
    for (i,gen) in pm_data["gen"]
        if gen["status_prev"] == 0 && gen["suqual"] == 1 && gen["gen_status"] == 0
            gen["gen_status"] = 1
            #println("$(i) / $(gen["source_id"]): $(0) -> $(1)")
        end
    end
    #println([gen["gen_status"] for (i,gen) in pm_data["gen"]])


    result = run_c2_opf_uc(pm_data, ACPPowerModel, nlp_solver, relax_integrality=true)

    if (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == OPTIMAL)

        committed = 0
        commitment_costs = 0.0
        improvement_uc = 0.0

        for (i,sol_gen) in result["solution"]["gen"]
            #println("$(sol_gen["gen_su"]), $(sol_gen["gen_sd"])")
            gen = pm_data["gen"][i]

            gen["gen_status"] = gen["status_prev"]

            if (sol_gen["gen_su"] >= 0.001 || sol_gen["gen_sd"] >= 0.001)
                println("$(i), $(sol_gen["gen_su"]), $(sol_gen["gen_sd"])")

                if sol_gen["gen_su"] > 0.0
                    if sol_gen["gen_su"] >= 0.90
                        println("  gen $(i) fixed to on, $(gen["sucost"])")
                        commitment_costs += gen["sucost"]
                        gen["gen_status"] = 1
                        committed += 1
                    end
                end

                if sol_gen["gen_sd"] > 0.0
                    if sol_gen["gen_sd"] >= 0.90
                        println("  gen $(i) fixed to off, $(gen["sdcost"])")
                        commitment_costs += gen["sdcost"]
                        gen["gen_status"] = 0
                        committed += 1
                    end
                end
            end
        end

        if committed > 0
            result = run_c2_opf_soft(pm_data, ACPPowerModel, nlp_solver)

            if (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == OPTIMAL)
                println("UC AC-OPF objective: $(result["objective"])")

                if result["objective"] > objective_best
                    improvement_uc = result["objective"] - objective_best
                    update_data!(pm_data, result["solution"])

                    correct_c2_solution!(pm_data)
                    write_c2_solution(pm_data, output_dir=output_dir)
                    objective_best = result["objective"]
                end
            else
                warn(_LOGGER, "run_c2_opf_soft failed after relaxed c2_opf_uc solve with status $(result["termination_status"])")
                improvement_uc = -1.0
                for (i,gen) in pm_data["gen"]
                    gen["gen_status"] = gen["status_prev"]
                end
            end
        else
            println("no commitment actions found")
        end

        println("COMMITTED, $(pm_data["name"]), $(committed), $(commitment_costs), $(improvement_uc)")
    else
        println("run_c2_opf_uc failed with $(result["termination_status"])")
    end


    ###### Prepare OTS Solution ######
    if scoring_method == 3 || scoring_method == 4
        println("")
        println("*** Start OTS Solve:")

        costs_updated = 0
        for (i,branch) in pm_data["branch"]
            if branch["status_prev"] == 0 && branch["swqual"] == 1 && branch["br_status"] == 0
                branch["br_status"] = 1
            end

            if branch["csw"] < 100.0
                branch["csw"] = 100.0
                costs_updated += 1
            end
        end

        if costs_updated > 0
            warn(_LOGGER, "updated $(costs_updated) branch costs to 100.0")
        end

        result = run_c2_ots_soft(pm_data, ACPPowerModel, nlp_solver, relax_integrality=true)

        if (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == OPTIMAL)
            println("OTS objective: $(result["objective"])")
            #print_summary(result["solution"])
            #return

            improvement_ots = 0.0
            switched = 0
            switch_costs = 0.0
            for (i,sol_branch) in result["solution"]["branch"]
                #println("$(sol_gen["gen_su"]), $(sol_gen["gen_sd"])")
                branch = pm_data["branch"][i]

                branch["br_status"] = branch["status_prev"]

                if sol_branch["z_branch_delta_abs"] >= 0.99
                    println("$(i), $(sol_branch["z_branch_delta_abs"]), $(branch["csw"])")
                    switch_costs += branch["csw"]
                    switched += 1

                    if sol_branch["br_status"] >= 0.5
                        println("  branch $(i) switched on")
                        branch["br_status"] = 1
                    else
                        println("  branch $(i) switched off")
                        branch["br_status"] = 0
                    end
                end
            end

            if switched > 0
                result = run_c2_opf_soft(pm_data, ACPPowerModel, nlp_solver)

                if (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == OPTIMAL)
                    println("Switched AC-OPF objective: $(result["objective"])")
                    #print_summary(result["solution"])

                    if result["objective"] > objective_best
                        improvement_ots =result["objective"] - objective_best
                        update_data!(pm_data, result["solution"])

                        correct_c2_solution!(pm_data)
                        write_c2_solution(pm_data, output_dir=output_dir)
                        objective_best = result["objective"]
                    end
                else
                    warn(_LOGGER, "run_c2_opf_soft failed after relaxed run_c2_ots_soft solve with status $(result["termination_status"])")
                    improvement_ots = -1.0
                    for (i,branch) in pm_data["branch"]
                        branch["br_status"] = branch["status_prev"]
                    end
                end
            else
                println("no switching actions found")
            end

            println("SWITCHING, $(pm_data["name"]), $(switched), $(switch_costs), $(improvement_ots)")
        else
            warn(_LOGGER, "run_c2_ots_soft failed with status $(result["termination_status"])")
        end
    end

    time_solve_basecase = time() - time_basecase_start
    @info("basecase eval time: $(time_solve_basecase)")

    println()
end



function code2(con_file::String, json_file::String, raw_file::String, reserved::String, time_limit::Int64, scoring_method::Int64, case_id::String, scenario_id::String; output_dir=".")
    ###### Load Data ######
    time_load_start = time()
    b0 = Base.gc_bytes()

    goc_data = parse_c2_files(raw_file, con_file, json_file, case_id=case_id, scenario_id=scenario_id)

    pm_data = build_c2_pm_model(goc_data)
    #print_summary(pm_data)

    data_bytes = Base.gc_bytes() - b0
    time_load = time() - time_load_start


    ###### Prepare Solution 2 ######
    time_contingencies_start = time()

    time_remaining = time_limit - (time() - time_load_start) - 20
    @info("contingency eval time limit: $(trunc(Int, time_remaining))")

    check_available_memory(data_bytes)

    contingencies = contingency_order(pm_data)

    workers = Distributed.workers()

    process_data = []

    cont_per_proc = length(contingencies)/length(workers)

    for p in 1:length(workers)
        cont_start = trunc(Int, ceil(1+(p-1)*cont_per_proc))
        cont_end = min(length(contingencies), trunc(Int,ceil(p*cont_per_proc)))
        pd = (
            pid = p,
            processes = length(workers),
            raw_file = raw_file,
            con_file = con_file,
            json_file = json_file,
            case_id = case_id,
            scenario_id = scenario_id,
            output_dir = output_dir,
            cont_range = cont_start:cont_end,
            time_limit = time_remaining
        )
        #println(pd)
        push!(process_data, pd)
    end

    @info("global time limit of $(time_remaining) seconds")
    for (i,pd) in enumerate(process_data)
        @info("worker task $(pd.pid): $(length(pd.cont_range)) / $(pd.cont_range)")
    end

    #for (i,pd) in enumerate(process_data)
    #    println(pd.pid)
    #    for cont in cont_order[pd.cont_range]
    #        println("$(cont.label)")
    #    end
    #end

    pmap(contingency_solver, process_data)

    time_contingencies = time() - time_contingencies_start
    @info("contingency eval time: $(time_contingencies)")

    println("CONT_DATA, $(pm_data["name"]), $(length(contingencies)), $(time_contingencies)")
    println()
end


function evaluation_summary(basepath::AbstractString, workers::Int, code1_runtime::Real, code2_runtime::Real)
    println("skip evaluation summary")
end


function contingency_solver(process_data)
    if length(process_data.cont_range) <= 0
        return
    end

    time_start = time()
    time_limit = time() + process_data.time_limit - 600.0
    #println(time_limit - time())

    PowerModels.silence()
    goc_data = parse_c2_files(process_data.raw_file, process_data.con_file, process_data.json_file, case_id=process_data.case_id, scenario_id=process_data.scenario_id)
    pm_data = build_c2_pm_model(goc_data)

    pm_sol = read_c2_solution(pm_data, output_dir=process_data.output_dir)
    PowerModels.update_data!(pm_data, pm_sol)
    update_previous!(pm_data)

    set_va_start_values!(pm_data)

    pm_data_base = deepcopy(pm_data)
    time_data = time() - time_start

    contingencies = contingency_order(pm_data)[process_data.cont_range]

    pm_data = deepcopy(pm_data_base)

    # write out basic solver
    for cont in contingencies
        if cont.type == "gen"
            println("working on: $(cont.label)")
            #pm_data = deepcopy(pm_data_base)
            pm_data["cont_label"] = cont.label

            basecase_gen = pm_data["gen"]["$(cont.idx)"]
            basecase_gen["present"] = false

            correct_c2_solution!(pm_data, contingency=true)
            write_c2_solution(pm_data, output_dir=process_data.output_dir, label=cont.label)

            basecase_gen["present"] = true

        elseif cont.type == "branch"
            println("working on: $(cont.label)")
            #pm_data = deepcopy(pm_data_base)
            pm_data["cont_label"] = cont.label

            basecase_branch = pm_data["branch"]["$(cont.idx)"]
            basecase_branch["present"] = false

            correct_c2_solution!(pm_data, contingency=true)
            write_c2_solution(pm_data, output_dir=process_data.output_dir, label=cont.label)

            basecase_branch["present"] = true
        else
            @assert("contingency type $(cont.type) not known")
        end
    end



    println("*****")
    println("Start CONT. OPT.")
    println("*****")

    contingencies_ordered = []

    pm_data_base_branch_flow = deepcopy(pm_data_base)
    ac_flows = calc_branch_flow_ac(pm_data_base_branch_flow)
    update_data!(pm_data_base_branch_flow, ac_flows)

    for cont in contingencies
        impact = 0.0

        if cont.type == "gen"
            gen_data = pm_data_base_branch_flow["gen"]["$(cont.idx)"]

            if gen_data["gen_status"] != 0
                impact = abs(gen_data["pg"] + gen_data["qg"]im)
            end
        elseif cont.type == "branch"
            branch_data = pm_data_base_branch_flow["branch"]["$(cont.idx)"]

            if branch_data["br_status"] != 0
                impact = abs(branch_data["pf"] + branch_data["qf"]im)
            end
        else
            @assert("contingency type $(cont.type) not known")
        end

        cont_data = (cont..., impact=impact)
        push!(contingencies_ordered, cont_data)
    end

    sort!(contingencies_ordered, by=(x) -> x.impact, rev=true)

    #println(contingencies_ordered)


    nlp_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

    solve_time_estimate = time_data
    for cont in contingencies_ordered
        if time() + 5.0*solve_time_estimate < time_limit
            time_solve_start = time()
            if cont.type == "gen"
                println("working on: $(cont.label)")
                pm_data["cont_label"] = cont.label

                basecase_gen = pm_data["gen"]["$(cont.idx)"]

                basecase_gen_status = basecase_gen["gen_status"]
                basecase_gen["gen_status"] = 0
                basecase_gen["present"] = false

                time_start = time()

                result = run_c2_opf_soft_ctg(pm_data, ACPPowerModel, nlp_solver)

                if (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
                    println("objective: $(result["objective"])")
                    #print_summary(result["solution"])

                    basecase_gen_pg = basecase_gen["pg"]
                    basecase_gen_qg = basecase_gen["qg"]

                    result["solution"]["gen"]["$(cont.idx)"] = Dict("pg" => 0.0, "qg" => 0.0)
                    update_data!(pm_data, result["solution"])

                    correct_c2_solution!(pm_data, contingency=true)
                    write_c2_solution(pm_data, output_dir=process_data.output_dir, label=cont.label)

                    basecase_gen["pg"] = basecase_gen_pg
                    basecase_gen["qg"] = basecase_gen_qg
                else
                    warn(_LOGGER, "run_c2_opf_soft_ctg failed with status $(result["termination_status"])")
                end

                basecase_gen["gen_status"] = basecase_gen_status
                basecase_gen["present"] = true

            elseif cont.type == "branch"
                println("working on: $(cont.label)")
                pm_data["cont_label"] = cont.label

                basecase_branch = pm_data["branch"]["$(cont.idx)"]

                basecase_branch_status = basecase_branch["br_status"]
                basecase_branch["br_status"] = 0
                basecase_branch["present"] = false

                time_start = time()

                result = run_c2_opf_soft_ctg(pm_data, ACPPowerModel, nlp_solver)

                if (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
                    println("objective: $(result["objective"])")
                    #print_summary(result["solution"])

                    update_data!(pm_data, result["solution"])

                    correct_c2_solution!(pm_data, contingency=true)
                    write_c2_solution(pm_data, output_dir=process_data.output_dir, label=cont.label)
                else
                    warn(_LOGGER, "run_c2_opf_soft_ctg failed with status $(result["termination_status"])")
                end

                basecase_branch["br_status"] = basecase_branch_status
                basecase_branch["present"] = true
            else
                @assert("contingency type $(cont.type) not known")
            end
            solve_time_estimate = time() - time_solve_start
        else
            warn(_LOGGER, "stopping contingency_solver due to time limit")
            break
        end
    end

end

