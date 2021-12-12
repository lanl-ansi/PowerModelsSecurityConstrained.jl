
function compute_c1_solution2_fast(con_file::String, inl_file::String, raw_file::String, rop_file::String, time_limit::Int, scoring_method::Int, network_model::String; output_dir::String="", scenario_id::String="none")
    time_data_start = time()
    goc_data = parse_c1_files(con_file, inl_file, raw_file, rop_file, scenario_id=scenario_id)
    network = build_c1_pm_model(goc_data)
    load_time = time() - time_data_start

    ###### Prepare Solution 2 ######

    time_contingencies_start = time()

    gen_cont_total = length(network["gen_contingencies"])
    branch_cont_total = length(network["branch_contingencies"])
    cont_total = gen_cont_total + branch_cont_total

    cont_order = contingency_order(network)

    #for cont in cont_order
    #    println(cont.label)
    #end

    workers = Distributed.workers()

    process_data = []

    cont_per_proc = cont_total/length(workers)

    for p in 1:length(workers)
        cont_start = trunc(Int, ceil(1+(p-1)*cont_per_proc))
        cont_end = min(cont_total, trunc(Int,ceil(p*cont_per_proc)))
        pd = (
            pid = p,
            processes = length(workers),
            con_file = con_file,
            inl_file = inl_file,
            raw_file = raw_file,
            rop_file = rop_file,
            scenario_id = scenario_id,
            output_dir = output_dir,
            cont_range = cont_start:cont_end,
        )
        #println(pd)
        push!(process_data, pd)
    end

    for (i,pd) in enumerate(process_data)
        info(LOGGER, "worker task $(pd.pid): $(length(pd.cont_range)) / $(pd.cont_range)")
    end

    #for (i,pd) in enumerate(process_data)
    #    println(pd.pid)
    #    for cont in cont_order[pd.cont_range]
    #        println("$(cont.label)")
    #    end
    #end


    solution2_files = pmap(c1_solution2_solver_fast, process_data)

    sort!(solution2_files)

    #println("pmap result: $(solution2_files)")

    time_contingencies = time() - time_contingencies_start
    info(LOGGER, "contingency eval time: $(time_contingencies)")

    info(LOGGER, "combine $(length(solution2_files)) solution2 files")
    c1_combine_files(solution2_files, "solution2.txt"; output_dir=output_dir)
    remove_c1_files(solution2_files)

    println("")

    data = [
        "----",
        "scenario id",
        "bus",
        "branch",
        "gen_cont",
        "branch_cont",
        "runtime (sec.)",
    ]
    println(join(data, ", "))

    data = [
        "DATA_SSS",
        goc_data.scenario,
        length(network["bus"]),
        length(network["branch"]),
        length(network["gen_contingencies"]),
        length(network["branch_contingencies"]),
        time_contingencies,
    ]
    println(join(data, ", "))

    write_c1_file_paths(goc_data.files; output_dir=output_dir)

    #println("")
    #write_evaluation_summary(goc_data, network, objective_lb=-Inf, load_time=load_time, contingency_time=time_contingencies, output_dir=output_dir)
end


@everywhere function c1_solution2_solver_fast(process_data)
    #println(process_data)
    time_data_start = time()
    PowerModels.silence()
    goc_data = parse_c1_files(
        process_data.con_file, process_data.inl_file, process_data.raw_file,
        process_data.rop_file, scenario_id=process_data.scenario_id)
    network = build_c1_pm_model(goc_data)

    sol = read_c1_solution1(network, output_dir=process_data.output_dir)
    PowerModels.update_data!(network, sol)
    time_data = time() - time_data_start

    for (i,bus) in network["bus"]
        if haskey(bus, "evhi")
            bus["vmax"] = bus["evhi"]
        end
        if haskey(bus, "evlo")
            bus["vmin"] = bus["evlo"]
        end
    end

    contingencies = contingency_order(network)[process_data.cont_range]

    pad_size = trunc(Int, ceil(log(10,process_data.processes)))
    padded_pid = lpad(string(process_data.pid), pad_size, "0")
    solution_filename = "solution2-$(padded_pid).txt"

    if length(process_data.output_dir) > 0
        solution_path = joinpath(process_data.output_dir, solution_filename)
    else
        solution_path = solution_filename
    end
    if isfile(solution_path)
        warn(LOGGER, "removing existing solution2 file $(solution_path)")
        rm(solution_path)
    end
    open(solution_path, "w") do sol_file
        # creates an empty file in the case of workers without contingencies
    end



    for cont in contingencies
        if cont.type == "gen"
            info(LOGGER, "working on: $(cont.label)")
            time_start = time()
            sol_cont = deepcopy(sol)
            debug(LOGGER, "contingency copy time: $(time() - time_start)")

            sol_cont["label"] = cont.label
            sol_cont["feasible"] = true
            sol_cont["cont_type"] = "gen"
            sol_cont["cont_comp_id"] = cont.idx

            sol_cont["delta"] = 0.0
            sol_cont["gen"]["$(cont.idx)"]["pg"] = 0.0
            sol_cont["gen"]["$(cont.idx)"]["qg"] = 0.0

            correct_c1_contingency_solution!(network, sol_cont)
            open(solution_path, "a") do sol_file
                sol2 = write_c1_solution2_contingency(sol_file, network, sol_cont)
            end
        elseif cont.type == "branch"
            info(LOGGER, "working on: $(cont.label)")
            time_start = time()
            sol_cont = deepcopy(sol)
            debug(LOGGER, "contingency copy time: $(time() - time_start)")
            sol_cont = deepcopy(sol)

            sol_cont["label"] = cont.label
            sol_cont["feasible"] = true
            sol_cont["cont_type"] = "branch"
            sol_cont["cont_comp_id"] = cont.idx

            sol_cont["delta"] = 0.0

            correct_c1_contingency_solution!(network, sol_cont)
            open(solution_path, "a") do sol_file
                sol2 = write_c1_solution2_contingency(sol_file, network, sol_cont)
            end
        else
            @assert("contingency type $(cont.type) not known")
        end
    end

    return solution_path
end
