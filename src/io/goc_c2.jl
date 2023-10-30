##### GOC C2 Data Tools #####

# controls tolerance to update hard constraints when writing solutions
const C2_FEAS_TOL = 1e-6

# controls tolerance of generator costs and load values that will be considered (slack value is 10000)
const C2_GEN_COST_TOL = 9999.0
const C2_LOAD_VAL_TOL = 9999.0


"""
Given an a directory path looks for and parses the files from ARPA-e GOC
Challenge 2 data format.
* `.raw`  network model
* `.con`  contingency set
* `.json` additional model data
"""
function parse_c2_case(case_directory)
    files = find_c2_case_files(case_directory)

    return parse_c2_files(files.raw, files.con, files.json, case_id=files.case_id, scenario_id=files.scenario_id)
end

function find_c2_case_files(case_directory; case_id=basename(dirname(case_directory)), scenario_id=basename(case_directory))
    if !isdir(case_directory)
        error(_LOGGER, "directory $(case_directory) not found")
    end

    if isfile(joinpath(case_directory, "case.raw"))
        raw_file = joinpath(case_directory, "case.raw")
    else
        raw_file = joinpath(case_directory, "case.RAW")
    end
    con_file = joinpath(case_directory, "case.con")
    json_file = joinpath(case_directory, "case.json")

    return (raw=raw_file, con=con_file, json=json_file, case_id=case_id, scenario_id=scenario_id)
end

function parse_c2_files(raw_file, con_file, json_file; case_id="none", scenario_id="none")
    files = Dict(
        "raw" => raw_file,
        "con" => con_file,
        "json" => json_file
    )

    info(_LOGGER, "Case Names")
    info(_LOGGER, "      case: $(case_id)")
    info(_LOGGER, "  scenario: $(scenario_id)")

    info(_LOGGER, "Parsing Files")
    info(_LOGGER, "   raw: $(files["raw"])")
    info(_LOGGER, "   con: $(files["con"])")
    info(_LOGGER, "  json: $(files["json"])")

    info(_LOGGER, "skipping power models data warnings")
    pm_logger_level = getlevel(getlogger(PowerModels))
    setlevel!(getlogger(PowerModels), "error")
    network_model = _PM.parse_file(files["raw"], import_all=true)
    setlevel!(getlogger(PowerModels), pm_logger_level)

    contingencies = parse_con_file(files["con"])

    json_data = Dict{String,Any}()
    open(files["json"], "r") do io
        json_data = JSON.parse(io)
    end
    # TODO Validate JSON data?

    return (files=files, case=case_id, scenario=scenario_id, network=network_model, contingencies=contingencies, json=json_data)
end



##### Transform GOC Data in PM Data #####

"""
Transforms files from ARPA-e GOC Challenge 2 data format in to the PowerModels
data format.  This consists of taking the data from multiple data structures
and putting into a network data dictionary.
"""
function build_c2_pm_model(goc_data)
    scenario = goc_data.scenario
    network = deepcopy(goc_data.network)

    network["name"] = "$(goc_data.case)-$(goc_data.scenario)"

    ##### General Helpers #####
    gen_lookup = Dict(tuple(gen["source_id"][2], gen["source_id"][3]) => gen for (i,gen) in network["gen"])

    load_lookup = Dict(tuple(load["source_id"][2], load["source_id"][3]) => load for (i,load) in network["load"])

    branch_lookup = Dict()
    for (i,branch) in network["branch"]
        if !branch["transformer"]
            branch_id = tuple(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][4]))
        else
            branch_id = tuple(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][5]))
            @assert branch["source_id"][4] == 0
            @assert branch["source_id"][6] == 0
        end
        branch_lookup[branch_id] = branch
    end


    ##### Network Data Checks #####
    @assert network["per_unit"]
    mva_base = network["baseMVA"]


    ##### Merge JSON Data #####
    json = goc_data.json

    network["pcblocks"] = json["pcblocks"]
    network["qcblocks"] = json["qcblocks"]
    network["scblocks"] = json["scblocks"]

    for (k,v) in json["systemparameters"]
        #println("$(k), $(v)")
        network[k] = v
    end

    for json_load in json["loads"]
        load_id = tuple(json_load["bus"], json_load["id"])
        if haskey(load_lookup, load_id)
            load = load_lookup[load_id]
            for (k,v) in json_load
                if k != "bus" && k != "id"
                    @assert(!haskey(load, k))
                    load[k] = v
                end
            end
        else
            warn(_LOGGER, "unable to find load $(load_id) from JSON file in RAW data.")
        end
    end

    for json_gen in json["generators"]
        gen_id = tuple(json_gen["bus"], json_gen["id"])
        if haskey(gen_lookup, gen_id)
            gen = gen_lookup[gen_id]
            for (k,v) in json_gen
                if k != "bus" && k != "id"
                    @assert(!haskey(gen, k))
                    gen[k] = v
                end
            end
        else
            warn(_LOGGER, "unable to find gen $(gen_id) from JSON file in RAW data.")
        end
    end

    for json_line in json["lines"]
        branch_id = tuple(json_line["origbus"], json_line["destbus"], json_line["id"])
        if haskey(branch_lookup, branch_id)
            branch = branch_lookup[branch_id]
            @assert(!branch["transformer"])
            for (k,v) in json_line
                if k != "origbus" && k != "destbus" && k != "id"
                    #println(branch_id, " ", k, " ", v)
                    @assert(!haskey(branch, k))
                    branch[k] = v
                end
            end
        else
            warn(_LOGGER, "unable to find line $(branch_id) from JSON file in RAW data.")
        end
    end

    for json_xfer in json["transformers"]
        branch_id = tuple(json_xfer["origbus"], json_xfer["destbus"], json_xfer["id"])
        if haskey(branch_lookup, branch_id)
            branch = branch_lookup[branch_id]
            @assert(branch["transformer"])
            for (k,v) in json_xfer
                if k != "origbus" && k != "destbus" && k != "id"
                    #println(branch_id, " ", k, " ", v)
                    @assert(!haskey(branch, k))
                    branch[k] = v
                end
            end
            #println(branch)
        else
            warn(_LOGGER, "unable to find transformer $(branch_id) from JSON file in RAW data.")
        end
    end


    ##### Transform System Cost Data #####
    network["p_delta_cost"] = cost_blocks_to_point_list(network["pcblocks"], mva_base, "pmax")
    network["q_delta_cost"] = cost_blocks_to_point_list(network["qcblocks"], mva_base, "qmax")
    network["sm_cost"] = cost_blocks_to_point_list(network["scblocks"], 1.0, "tmax") # given in percentage

    network["power_vio_limit"] = 10.0
    #network["power_vio_limit"] = 10000000.0
    network["p_delta_cost_approx"] = calc_c2_linear_cost_approximation(network["p_delta_cost"], 10.0)
    network["q_delta_cost_approx"] = calc_c2_linear_cost_approximation(network["q_delta_cost"], 10.0)

    network["sm_vio_limit"] = 0.2
    network["sm_cost_approx"] = calc_c2_linear_cost_approximation(network["sm_cost"], 0.2)



    ##### Impedance Correction Updates #####
    ict = Dict{Int,Any}()
    network["impedance_correction_table"] = ict

    if haskey(network, "impedance correction")
        for (i,ict_raw) in network["impedance correction"]
            #println(ict_raw)

            ict_index = ict_raw["i"]

            points = []
            for i in 1:11
                t_key = "t$(i)"
                f_key = "f$(i)"
                if haskey(ict_raw, t_key) && haskey(ict_raw, f_key)
                    point = (t=ict_raw[t_key], f=ict_raw[f_key])

                    if !(isapprox(point.t, 0.0) && isapprox(point.f, 0.0))
                        push!(points, point)
                    else
                        #warn(_LOGGER, "skipping point $(point)")
                    end
                else
                    break
                end
            end

            if length(points) <= 1
                error(_LOGGER, "no valid points found for impedance correction data id $(ict_index)")
            end
            if length(points) == 2
                error(_LOGGER, "only one valid points found for impedance correction data id $(ict_index), at least two are required")
            end

            if haskey(ict, ict_index)
                error(_LOGGER, "Impedance Correction table id $(ict_index) defined twice")
            end

            ict[ict_index] = points
        end

        delete!(network, "impedance correction")
    end

    ##### Bus Data Updates #####
    for (i,bus) in network["bus"]
        bus["present"] = true
    end


    ##### Generator Data Updates #####
    for (i,gen) in network["gen"]
        if !haskey(gen, "cblocks") || length(gen["cblocks"]) < 1
            warn(_LOGGER, "cost model data missing from gen $(i), $(gen["source_id"])")
            cblock = Dict{String,Any}(
                "pmax" => gen["pmax"],
                "c" => 0.0
            )
            push!(gen["cblocks"], cblock)
        end

        if length(gen["cblocks"]) > 1 && any(cblock["c"] >= C2_GEN_COST_TOL for cblock in gen["cblocks"])
            warn(_LOGGER, "update cost model data on gen $(i), $(gen["source_id"]), to remove a cost block that is above $(C2_GEN_COST_TOL)")
            gen["cblocks"] = [cblock for cblock in gen["cblocks"] if cblock["c"] <= C2_GEN_COST_TOL]
            cost_pmax = sum(cblock["pmax"] for cblock in gen["cblocks"])
            gen["pmax"] = min(gen["pmax"], cost_pmax)
        end

        pmax = gen["pmax"]*mva_base
        for cblock in gen["cblocks"]
             if cblock["pmax"] > pmax && cblock["pmax"] < pmax*1.1 
                warn(_LOGGER, "update cost model data on gen $(i), $(gen["source_id"]), reduce cost block max from $(cblock["pmax"]) to $(pmax*1.1)")
                if isapprox(pmax, 0.0)
                    cblock["pmax"] = 1.0
                else
                    cblock["pmax"] = pmax*1.1
                end
            end
        end

        point_list = cost_blocks_to_point_list(gen["cblocks"], mva_base, "pmax")

        gen["model"] = 1
        gen["ncost"] = div(length(point_list), 2)
        gen["cost"] = point_list

        gen["status_prev"] = gen["gen_status"]
        gen["pg_prev"] = gen["pg"]
        gen["qg_prev"] = gen["qg"]

        gen["prumax"] = gen["prumax"]/mva_base
        gen["prdmax"] = gen["prdmax"]/mva_base
        gen["prumaxctg"] = gen["prumaxctg"]/mva_base
        gen["prdmaxctg"] = gen["prdmaxctg"]/mva_base
        #println(i, " ", point_list)

        gen["present"] = true
    end


    ##### Load Data Updates #####
    for (i,load) in network["load"]
        if !haskey(load, "cblocks") || length(load["cblocks"]) < 1
            warn(_LOGGER, "cost model data missing from load $(i), $(load["source_id"])")
            cblock = Dict{String,Any}(
                "pmax" => load["pd"]*load["tmax"],
                "c" => 0.0
            )
            push!(load["cblocks"], cblock)
        end

        if length(load["cblocks"]) > 1 && any(cblock["c"] >= C2_LOAD_VAL_TOL for cblock in load["cblocks"])
            warn(_LOGGER, "value model data on load $(i), $(load["source_id"]) is above $(C2_LOAD_VAL_TOL)")
        end

        pmax = load["pd"]*load["tmax"]*mva_base
        for cblock in load["cblocks"]
            if cblock["pmax"] > pmax && pmax*1.1 < cblock["pmax"]
                warn(_LOGGER, "update cost model data on load $(i), $(load["source_id"]), reduce cost block max from $(cblock["pmax"]) to $(pmax*1.1)")
                cblock["pmax"] = pmax*1.1
                if isapprox(pmax, 0.0)
                    cblock["pmax"] = 1.0
                else
                    cblock["pmax"] = pmax*1.1
                end
            end
        end

        point_list = cost_blocks_to_point_list(load["cblocks"], mva_base, "pmax", decreasing=true)

        load["model"] = 1
        load["ncost"] = div(length(point_list), 2)
        load["cost"] = point_list

        load["pd_nominal"] = load["pd"]
        load["qd_nominal"] = load["qd"]
        load["pd_max"] = load["pd_nominal"]*load["tmax"]
        load["pd_min"] = load["pd_nominal"]*load["tmin"]

        # needed for reading in t-values from solutions files
        load["pd_prev"] = load["pd"]
        load["qd_prev"] = load["qd"]

        load["prumax"] = load["prumax"]/mva_base
        load["prdmax"] = load["prdmax"]/mva_base
        load["prumaxctg"] = load["prumaxctg"]/mva_base
        load["prdmaxctg"] = load["prdmaxctg"]/mva_base

        load["present"] = (load["status"] != 0)
    end


    ##### Shunt Data Updates #####
    for (i,shunt) in network["shunt"]
        #println(shunt["source_id"])
        if shunt["source_id"][1] == "switched shunt"
            #@assert shunt["source_id"][3] == 0
            @assert shunt["gs"] == 0.0
            shunt["dispatchable"] = true

            blocks = Tuple{Int64,Float64}[]
            for (n_name,b_name) in [("n1","b1"),("n2","b2"),("n3","b3"),("n4","b4"),("n5","b5"),("n6","b6"),("n7","b7"),("n8","b8")]
                if !isapprox(shunt[b_name], 0.0)
                    push!(blocks, (shunt[n_name], shunt[b_name]/mva_base))
                end
            end
            shunt["blocks"] = blocks

            #println(shunt["source_id"], " ", blocks)
            if length(blocks) > 0
                b_value, combination = calc_c2_closest_discrete_combination(shunt["bs"], blocks)
                shunt["bs"] = b_value
                shunt["xst"] = combination
            else
                #warn(_LOGGER, "switching data missing from shunt $(i)")
                @assert shunt["bs"] == 0.0
                shunt["xst"] = []
            end

            bmin = 0.0
            bmax = 0.0
            for (num,sup) in blocks
                if sup <= 0.0
                    bmin += num*sup
                else
                    bmax += num*sup
                end
            end
            shunt["bmin"] = bmin
            shunt["bmax"] = bmax

            #println("$(blocks), $(shunt["bmin"]), $(shunt["bmax"])")

            shunt["present"] = (shunt["status"] != 0)
        else
            shunt["dispatchable"] = false
            shunt["present"] = true
        end
    end


    ##### Branch Data Updates #####
    for (i,branch) in network["branch"]
        branch["present"] = true
        branch["status_prev"] = branch["br_status"]

        if branch["transformer"]
            control_mode = branch["control_mode"] = branch["cod1"]

            branch["br_r_nominal"] = branch["br_r"]
            branch["br_x_nominal"] = branch["br_x"]

            if branch["source_id"][2] == branch["t_bus"] && control_mode != 0
                warn(_LOGGER, "reverting parallel transformer branch $(branch["index"]) to data orientation so that tap indexing will be correct")
                branch_orginal = copy(branch)
                branch["f_bus"] = branch_orginal["t_bus"]
                branch["t_bus"] = branch_orginal["f_bus"]
                branch["g_to"] = branch_orginal["g_fr"] .* branch_orginal["tap"]'.^2
                branch["b_to"] = branch_orginal["b_fr"] .* branch_orginal["tap"]'.^2
                branch["g_fr"] = branch_orginal["g_to"] ./ branch_orginal["tap"]'.^2
                branch["b_fr"] = branch_orginal["b_to"] ./ branch_orginal["tap"]'.^2
                branch["tap"] = 1 ./ branch_orginal["tap"]
                branch["br_r"] = branch_orginal["br_r"] .* branch_orginal["tap"]'.^2
                branch["br_x"] = branch_orginal["br_x"] .* branch_orginal["tap"]'.^2
                branch["shift"] = -branch_orginal["shift"]
                branch["angmin"] = -branch_orginal["angmax"]
                branch["angmax"] = -branch_orginal["angmin"]
            end

            if control_mode == 0
                # base case do nothing
            elseif control_mode == 1 || control_mode == -1
                branch["tm_max"] = branch["rma1"]
                branch["tm_min"] = branch["rmi1"]
                #println(branch["source_id"], " ", branch["tap"], " ", branch["rmi1"], " ", branch["rma1"])
                branch["tm_steps"] = trunc(Int, branch["ntp1"])

                if branch["tap"] < branch["tm_min"]
                    warn(_LOGGER, "transformer branch $(i) tap updated to be in bounds, $(branch["tap"]) -> $(branch["tm_min"])")
                    branch["tap"] = branch["tm_min"]
                end
                if branch["tap"] > branch["tm_max"]
                    warn(_LOGGER, "transformer branch $(i) tap updated to be in bounds, $(branch["tap"]) -> $(branch["tm_max"])")
                    branch["tap"] = branch["tm_max"]
                end

                tm_step, tm_value = closest_discrete_value(branch["tap"], branch["tm_min"], branch["tm_max"], branch["tm_steps"])
                branch["tm_step"] = tm_step
                branch["tap"] = tm_value

                if branch["tab1"] != 0
                    ict = network["impedance_correction_table"][branch["tab1"]]
                    correct_transformer_impedance!(branch, ict, branch["tap"])
                end

                #println(branch["source_id"], " ", branch["tap"], " ", branch["tm_step"], " ", value)
            elseif control_mode == 3 || control_mode == -3
                branch["ta_max"] = deg2rad(branch["rma1"])
                branch["ta_min"] = deg2rad(branch["rmi1"])
                branch["ta_steps"] = trunc(Int, branch["ntp1"])

                if branch["shift"] < branch["ta_min"]
                    warn(_LOGGER, "transformer branch $(i) shift updated to be in bounds, $(branch["shift"]) -> $(branch["ta_min"])")
                    branch["shift"] = branch["ta_min"]
                end
                if branch["shift"] > branch["ta_max"]
                    warn(_LOGGER, "transformer branch $(i) shift updated to be in bounds, $(branch["shift"]) -> $(branch["ta_max"])")
                    branch["shift"] = branch["ta_max"]
                end

                ta_step, ta_value = closest_discrete_value(branch["shift"], branch["ta_min"], branch["ta_max"], branch["ta_steps"])
                branch["ta_step"] = ta_step
                branch["shift"] = ta_value

                if branch["tab1"] != 0
                    ict = network["impedance_correction_table"][branch["tab1"]]
                    correct_transformer_impedance!(branch, ict, rad2deg(branch["shift"]))
                end

                #println("$(branch["source_id"]) $(branch["ta_step"]) $(rad2deg(branch["shift"]))")
            else
                # same as base case
                warn(_LOGGER, "transformer branch $(i) given invalid control mode value of $(control_mode), only values of -3,-1,0,1,3 are supported")
            end
        end
    end


    ##### Add Contingency Lists #####
    generator_ids = []
    branch_ids = []

    for (i,cont) in enumerate(goc_data.contingencies)
        if cont["component"] == "branch"
            branch_id = (cont["i"], cont["j"], cont["ckt"])
            pm_branch = branch_lookup[branch_id]
            push!(branch_ids, (idx=pm_branch["index"], label=cont["label"], type="branch"))

        elseif cont["component"] == "generator"
            gen_id = (cont["i"], cont["id"])
            pm_gen = gen_lookup[gen_id]
            push!(generator_ids, (idx=pm_gen["index"], label=cont["label"], type="gen"))

        else
            error(_LOGGER, "unrecognized contingency component type $(cont["component"]) at contingency $(i)")
        end
    end

    network["branch_contingencies"] = branch_ids
    network["gen_contingencies"] = generator_ids

    network["branch_contingencies_active"] = []
    network["gen_contingencies_active"] = []



    ##### Fix Broken Data #####
    _PM.correct_cost_functions!(network)

    return network
end


function cost_blocks_to_point_list(cost_blocks, mva_base, domain_key; fixed_cost=0.0, decreasing=false)
    cummulative_power = 0.0
    cummulative_cost = fixed_cost

    point_list = Float64[]
    push!(point_list, cummulative_power)
    push!(point_list, cummulative_cost)

    cost_blocks = [b for b in cost_blocks]
    if !decreasing
        sort!(cost_blocks, by=(x) -> x["c"])
    else
        sort!(cost_blocks, by=(x) -> x["c"], rev=true)
    end
    for block in cost_blocks
        cummulative_power += block[domain_key]/mva_base
        cummulative_cost += block["c"]*block[domain_key]
        push!(point_list, cummulative_power)
        push!(point_list, cummulative_cost)
    end

    return point_list
end

"computes a pesimistic linear approximation from pwl cost points"
function calc_c2_linear_cost_approximation(point_list, target_value)
    active_segement = 0
    for i in 0:(div(length(point_list), 2) - 1)
        active_segement = i
        if target_value >= point_list[2*i+1] && target_value <= point_list[2*i+3]
            break
        end
    end
    #println(active_segement)
    x1 = point_list[2*active_segement+1]
    y1 = point_list[2*active_segement+2]
    x2 = point_list[2*active_segement+3]
    y2 = point_list[2*active_segement+4]

    #println(x1, " ", y1, " - ", x2, " ", y2)
    slope = (y1 - y2)/(x1 - x2)
    intercept = y1 - slope*(x1)
    target_cost = slope*(target_value) + intercept

    #println(target_value, " ", target_cost)
    return round(target_cost/target_value)
end


"finds the discrete settings that are closest to a given floating point value"
function calc_c2_closest_discrete_combination(value::Real, blocks::Array{Tuple{Int64,Float64},1})
    for block in blocks
        @assert(block[1] >= 0)
    end

    combinations = []
    _c2_comp_combinations!(combinations, blocks, Int64[])

    num_blocks = length(blocks)

    closest_comb = [0 for i in 1:num_blocks]
    closest_value = 0.0

    for combination in combinations
        combination_value = sum(blocks[i][2]*combination[i] for i in 1:num_blocks)
        #return target_step, closest_value
        #println(value, " ", combination_value, " ", combination)
        if abs(combination_value - value) <= abs(closest_value - value)
            if isapprox(abs(combination_value - value), abs(closest_value - value)) && sum(combination) >= sum(closest_comb)
                continue
            end
            closest_comb = combination
            closest_value = combination_value
        end
    end

    return closest_value, closest_comb
end

function _c2_comp_combinations!(combinations, blocks::Array{Tuple{Int64,Float64},1}, partial_assignment)
    if length(blocks) < 1
        push!(combinations, partial_assignment)
    else
        for i in 0:blocks[1][1]
            npa = vcat(partial_assignment, i)
            _c2_comp_combinations!(combinations, blocks[2:end], npa)
        end
    end
end


"checks feasibility criteria of network solution, corrects when possible"
function correct_c2_solution!(network; contingency=false)

    delta_r = network["deltar"]
    if contingency
        delta_r = network["deltarctg"]
    end

    # default value is required for correctness of max and mean computations when no changes are made 
    vm_changes = [0.0]
    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            if bus["vm"] > bus["vmax"]
                warn(_LOGGER, "update vm on bus $(i) to be in bounds $(bus["vm"]) -> $(bus["vmax"])")
                push!(vm_changes, bus["vm"] - bus["vmax"])
                bus["vm"] = bus["vmax"]
            end

            if bus["vm"] < bus["vmin"]
                warn(_LOGGER, "update vm on bus $(i) to be in bounds $(bus["vm"]) -> $(bus["vmin"])")
                push!(vm_changes, bus["vmin"] - bus["vm"])
                bus["vm"] = bus["vmin"]
            end
        else
            bus["vm"] = 0.0
            bus["va"] = 0.0
        end
    end

    load_changes = [0.0]
    for (i,load) in network["load"]
        if load["status"] != 0
            prdmax = load["prdmax"]
            prumax = load["prdmax"]
            if contingency
                prdmax = load["prdmaxctg"]
                prumax = load["prumaxctg"]
            end

            if load["pd"] < load["pd_prev"] - prdmax*delta_r - C2_FEAS_TOL
                new_pd = load["pd_prev"] - prdmax*delta_r
                warn(_LOGGER, "update pd on load $(i) to be in ramping bounds $(load["pd"]) -> $(new_pd)")
                load["pd"] = new_pd
            end

            if load["pd"] > load["pd_prev"] + prumax*delta_r + C2_FEAS_TOL
                new_pd = load["pd_prev"] + prumax*delta_r
                warn(_LOGGER, "update pd on load $(i) to be in ramping bounds $(load["pd"]) -> $(new_pd)")
                load["pd"] = new_pd
            end

            if load["pd"] > load["pd_max"]
                warn(_LOGGER, "update power on load $(i) to be in bounds $(load["pd"]) -> $(load["pd_max"])")
                push!(load_changes, load["pd"] - load["pd_max"])
                load["pd"] = load["pd_max"]
            end
            if load["pd"] < load["pd_min"]
                warn(_LOGGER, "update power on load $(i) to be in bounds $(load["pd"]) -> $(load["pd_min"])")
                push!(load_changes, load["pd_min"] - load["pd"])
                load["pd"] = load["pd_min"]
            end

            if isapprox(load["pd_nominal"], 0.0)
                if load["qd_nominal"] >= 0
                    qd_max = load["qd_nominal"]*load["tmax"]
                    qd_min = load["qd_nominal"]*load["tmin"]
                else
                    qd_max = load["qd_nominal"]*load["tmin"]
                    qd_min = load["qd_nominal"]*load["tmax"]
                end

                if load["qd"] > qd_max
                    warn(_LOGGER, "update reactive power on load $(i) to be in bounds $(load["qd"]) -> $(qd_max) : $(qd_min) to $(qd_max)")
                    load["qd"] = qd_max
                end
                if load["qd"] < qd_min
                    warn(_LOGGER, "update reactive power on load $(i) to be in bounds $(load["qd"]) -> $(qd_min) : $(qd_min) to $(qd_max)")
                    load["qd"] = qd_min
                end
            end

        end
    end

    pg_changes = [0.0]
    qg_changes = [0.0]
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            prdmax = gen["prdmax"]
            prumax = gen["prdmax"]
            if contingency
                prdmax = gen["prdmaxctg"]
                prumax = gen["prumaxctg"]
            end

            pg_prev = gen["pg_prev"]
            if gen["status_prev"] == 0 # starting up
               pg_prev = gen["pmin"]
               #warn(_LOGGER, "generator $(i) startung up setting previous value to $(pg_prev)")
            end

            if gen["pg"] < pg_prev - prdmax*delta_r - C2_FEAS_TOL
                new_pg = pg_prev - prdmax*delta_r
                warn(_LOGGER, "update pg on gen $(i) to be in ramping bounds $(gen["pg"]) -> $(new_pg)")
                gen["pg"] = new_pg
            end

            if gen["pg"] > pg_prev + prumax*delta_r + C2_FEAS_TOL
                new_pg = pg_prev + prumax*delta_r
                warn(_LOGGER, "update pg on gen $(i) to be in ramping bounds $(gen["pg"]) -> $(new_pg)")
                gen["pg"] = new_pg
            end

            if gen["pg"] > gen["pmax"]
                warn(_LOGGER, "update pg on gen $(i) to be in bounds $(gen["pg"]) -> $(gen["pmax"])")
                push!(pg_changes, gen["pg"] - gen["pmax"])
                gen["pg"] = gen["pmax"]
            end
            if gen["pg"] < gen["pmin"]
                warn(_LOGGER, "update pg on gen $(i) to be in bounds $(gen["pg"]) -> $(gen["pmin"])")
                push!(pg_changes, gen["pmin"] - gen["pg"])
                gen["pg"] = gen["pmin"]
            end

            if gen["qg"] > gen["qmax"]
                warn(_LOGGER, "update qg on gen $(i) to be in bounds $(gen["qg"]) -> $(gen["qmax"])")
                push!(qg_changes, gen["qg"] - gen["qmax"])
                gen["qg"] = gen["qmax"]
            end
            if gen["qg"] < gen["qmin"]
                warn(_LOGGER, "update qg on gen $(i) to be in bounds $(gen["qg"]) -> $(gen["qmin"])")
                push!(qg_changes, gen["qmin"] - gen["qg"])
                gen["qg"] = gen["qmin"]
            end
        else
            gen["pg"] = 0.0
            gen["qg"] = 0.0
        end
    end

end


function write_c2_solution(network; output_dir="", label="BASECASE")
    solution_file = "solution_$(label).txt"
    if length(output_dir) > 0
        solution_path = joinpath(output_dir, solution_file)
    else
        solution_path = solution_file
    end

    open(solution_path, "w") do sol1

        write(sol1, "--bus section\n")
        write(sol1, "i, v, theta\n")
        for (i,bus) in network["bus"]
            if bus["present"]
                write(sol1, "$(bus["index"]), $(bus["vm"]), $(bus["va"])\n")
            end
        end

        write(sol1, "--load section\n")
        write(sol1, "i, id, t\n")
        for (i,load) in network["load"]
            if load["present"]
                bus_index = load["source_id"][2]
                load_id = load["source_id"][3]
                if !isapprox(load["pd_nominal"], 0.0)
                    load_t = load["pd"]/load["pd_nominal"]
                elseif !isapprox(load["qd_nominal"], 0.0)
                    load_t = load["qd"]/load["qd_nominal"]
                else
                    load_t = (load["tmax"] + load["tmin"])/2.0
                end
                write(sol1, "$(bus_index), $(load_id), $(load_t)\n")
            end
        end

        write(sol1, "--generator section\n")
        write(sol1, "i, id, p, q, x\n")
        for (i,gen) in network["gen"]
            if gen["present"]
                bus_index = gen["source_id"][2]
                gen_id = gen["source_id"][3]
                write(sol1, "$(bus_index), $(gen_id), $(gen["pg"]), $(gen["qg"]), $(Int(gen["gen_status"] != 0))\n")
            end
        end

        write(sol1, "--line section\n")
        write(sol1, "iorig, idest, id, x\n")
        for (i,branch) in network["branch"]
            if branch["present"] && !branch["transformer"]
                bus_fr_index = branch["source_id"][2]
                bus_to_index = branch["source_id"][3]
                branch_ckt = branch["source_id"][4]
                write(sol1, "$(bus_fr_index), $(bus_to_index), $(branch_ckt), $(Int(branch["br_status"] != 0))\n")
            end
        end

        write(sol1, "--transformer section\n")
        write(sol1, "iorig, idest, id, x, xst\n")
        for (i,branch) in network["branch"]
            if branch["present"] && branch["transformer"]
                bus_fr_index = branch["source_id"][2]
                bus_to_index = branch["source_id"][3]
                branch_ckt = branch["source_id"][5]
                step_setting = 0

                if branch["control_mode"] == 0
                    # base case do nothing
                elseif branch["control_mode"] == 1 || branch["control_mode"] == -1
                    step_setting = branch["tm_step"]
                elseif branch["control_mode"] == 3 || branch["control_mode"] == -3
                    step_setting = branch["ta_step"]
                else
                end

                write(sol1, "$(bus_fr_index), $(bus_to_index), $(branch_ckt), $(Int(branch["br_status"] != 0)), $(step_setting)\n")
            end
        end

        write(sol1, "--switched shunt section\n")
        write(sol1, "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8\n")
        for (i,shunt) in network["shunt"]
            if shunt["present"] && shunt["dispatchable"]
                bus_index = shunt["source_id"][2]
                items = ["$(bus_index)"]
                for b in shunt["xst"]
                    push!(items, "$(b)")
                end
                write(sol1, "$(join(items, ", "))\n")
            end
        end

    end

end


function remove_c2_solution_files(; output_dir=".")
    for file in readdir(output_dir)
        if startswith(file, "solution_") && endswith(file, ".txt")
            info(_LOGGER, "deleting solution file: $(file) in $(output_dir)")
            rm(joinpath(output_dir, file))
        end
    end
end


function read_c2_solution(network; output_dir="", label="BASECASE")
    solution_file = "solution_$(label).txt"
    if length(output_dir) > 0
        solution_path = joinpath(output_dir, solution_file)
    else
        solution_path = solution_file
    end

    info(_LOGGER, "loading solution file: $(solution_path)")
    goc_sol = parse_c2_solution_file(solution_path)

    info(_LOGGER, "converting GOC solution to PowerModels solution")
    pm_sol = build_c2_pm_solution(network, goc_sol)

    return pm_sol
end


function parse_c2_solution_file(file::String)
    open(file) do io
        return parse_c2_solution_file(io)
    end
end

function parse_c2_solution_file(io::IO)
    bus_data_list = []
    load_data_list = []
    generator_data_list = []
    line_data_list = []
    transformer_data_list = []
    switched_shunt_data_list = []

    lines = readlines(io)

    # skip bus list header section
    idx = 1

    separator_count = 0
    skip_next = false

    while idx <= length(lines)
        line = lines[idx]
        if length(strip(line)) == 0
            warn(_LOGGER, "skipping blank line in solution file ($(idx))")
        elseif skip_next
            skip_next = false
        elseif startswith(strip(line), "--")
            separator_count += 1
            skip_next = true
        else
            if separator_count == 1
                bus_data = _parse_comma_seperated_line(line,
                    (:i, :v, :theta),
                    [Int, Float64, Float64]
                )
                push!(bus_data_list, bus_data)
            elseif separator_count == 2
                load_data = _parse_comma_seperated_line(line,
                    (:i, :id, :t),
                    [Int, String, Float64]
                )
                push!(load_data_list, load_data)
            elseif separator_count == 3
                generator_data = _parse_comma_seperated_line(line,
                    (:i, :id, :p, :q, :x),
                    [Int, String, Float64, Float64, Int]
                )
                push!(generator_data_list, generator_data)
            elseif separator_count == 4
                line_data = _parse_comma_seperated_line(line,
                    (:iorig, :idest, :id, :x),
                    [Int, Int, String, Int]
                )
                push!(line_data_list, line_data)
            elseif separator_count == 5
                transformer_data = _parse_comma_seperated_line(line,
                    (:iorig, :idest, :id, :x, :xst),
                    [Int, Int, String, Int, Int]
                )
                push!(transformer_data_list, transformer_data)
            elseif separator_count == 6
                parts = split(line, ",")
                @assert length(parts) >= 1

                parts = [parse(Int, v) for v in parts]
                switched_shunt_data = (i=parts[1], xst=[v for v in parts[2:end]])
                push!(switched_shunt_data_list, switched_shunt_data)
            else
                warn(_LOGGER, "skipping line in solution file ($(idx)): $(line)")
            end
        end
        idx += 1
    end

    return (bus=bus_data_list, load=load_data_list,
        generator=generator_data_list, line=line_data_list,
        transformer=transformer_data_list,
        switched_shunt=switched_shunt_data_list)
end

function Base.parse(type::Type{String}, value::AbstractString)
    return strip(value)
end

function _parse_comma_seperated_line(line, names, types)
    @assert length(names) >= length(types)

    parts = split(line, ",")
    @assert length(parts) >= length(names)

    data = NamedTuple{names}([parse(types[i], parts[i]) for i in 1:length(names)])

    return data
end


function build_c2_pm_solution(network, goc_sol)
    bus_lookup = Dict(parse(Int, bus["source_id"][2]) => bus for (i,bus) in network["bus"])
    load_lookup = Dict((load["source_id"][2], strip(load["source_id"][3])) => load for (i,load) in network["load"])
    gen_lookup = Dict((gen["source_id"][2], strip(gen["source_id"][3])) => gen for (i,gen) in network["gen"])

    branch_lookup = Dict()
    for (i,branch) in network["branch"]
        if !branch["transformer"]
            branch_lookup[(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][4]))] = branch
        else
            branch_lookup[(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][5]))] = branch
        end
    end

    shunt_lookup = Dict{Int,Any}()
    for (i,shunt) in network["shunt"]
        if shunt["source_id"][1] == "switched shunt"
            #@assert shunt["source_id"][3] == 0
            shunt_lookup[shunt["source_id"][2]] = shunt
        end
    end


    mva_base = network["baseMVA"]

    bus_data = Dict{String,Any}()
    for bus_sol in goc_sol.bus
        pm_bus = bus_lookup[bus_sol.i]
        bus_data["$(pm_bus["index"])"] = Dict(
            "vm" => bus_sol.v,
            "va" => bus_sol.theta
        )
    end

    load_data = Dict{String,Any}()
    for load_sol in goc_sol.load
        pm_load = load_lookup[(load_sol.i, load_sol.id)]
        load_data["$(pm_load["index"])"] = Dict(
            "pd" => pm_load["pd_nominal"]*load_sol.t,
            "qd" => pm_load["qd_nominal"]*load_sol.t
        )
    end

    gen_data = Dict{String,Any}()
    for gen_sol in goc_sol.generator
        pm_gen = gen_lookup[(gen_sol.i, gen_sol.id)]
        gen_data["$(pm_gen["index"])"] = Dict(
            "pg" => gen_sol.p,
            "qg" => gen_sol.q,
            "gen_status" => gen_sol.x
        )
    end

    shunt_data = Dict{String,Any}()
    for shunt_sol in goc_sol.switched_shunt
        pm_shunt = shunt_lookup[shunt_sol.i]

        bs_value = 0.0
        for (i,v) in enumerate(shunt_sol.xst)
            b_key = "b$(i)"
            if haskey(pm_shunt, b_key)
                bs_value += pm_shunt[b_key]/mva_base*v
            else
                warn(_LOGGER, "solution file had swtiched shunt parameters that were not present in the data file!")
            end
        end

        shunt_data["$(pm_shunt["index"])"] = Dict(
            "xst" => shunt_sol.xst,
            "bs" => bs_value
        )
    end

    branch_data = Dict{String,Any}()
    for line_sol in goc_sol.line
        pm_branch = branch_lookup[(line_sol.iorig, line_sol.idest, line_sol.id)]
        branch_data["$(pm_branch["index"])"] = Dict(
            "br_status" => line_sol.x
        )
    end
    for transformer_sol in goc_sol.transformer
        pm_branch = branch_lookup[(transformer_sol.iorig, transformer_sol.idest, transformer_sol.id)]

        #println(transformer_sol, " ", pm_branch["control_mode"])

        xfer_data = Dict{String,Any}(
            "br_status" => transformer_sol.x,
        )

        if pm_branch["control_mode"] == 0
            # base case do nothing
        elseif pm_branch["control_mode"] == 1 || pm_branch["control_mode"] == -1
            xfer_data["tm_step"] = transformer_sol.xst
            xfer_data["tap"] = tm_step_to_value(pm_branch, transformer_sol.xst)
            if pm_branch["tab1"] != 0
                ict = network["impedance_correction_table"][pm_branch["tab1"]]
                settings = comp_transformer_impedance!(pm_branch, ict, pm_branch["tap"])
                xfer_data["br_r"] = settings["br_r"]
                xfer_data["br_x"] = settings["br_x"]
                #println(xfer_data)
            end
        elseif pm_branch["control_mode"] == 3 || pm_branch["control_mode"] == -3
            xfer_data["ta_step"] = transformer_sol.xst
            xfer_data["shift"] = ta_step_to_value(pm_branch, transformer_sol.xst)
            if pm_branch["tab1"] != 0
                ict = network["impedance_correction_table"][pm_branch["tab1"]]
                settings = comp_transformer_impedance!(pm_branch, ict, rad2deg(pm_branch["shift"]))
                xfer_data["br_r"] = settings["br_r"]
                xfer_data["br_x"] = settings["br_x"]
                #println(xfer_data)
            end
        else
        end

        branch_data["$(pm_branch["index"])"] = xfer_data
    end

    solution = Dict(
        "per_unit" => true,
        "bus" => bus_data,
        "load" => load_data,
        "shunt" => shunt_data,
        "gen" => gen_data,
        "branch" => branch_data
    )

    return solution
end

function update_previous!(data)
    for (i,load) in data["load"]
        load["pd_prev"] = load["pd"]
        load["qd_prev"] = load["qd"]
    end

    for (i,gen) in data["gen"]
        gen["status_prev"] = gen["gen_status"]
        gen["pg_prev"] = gen["pg"]
        gen["qg_prev"] = gen["qg"]
    end
end

function tm_step_to_value(branch, tm_step::Int)
    mid_point = (branch["tm_max"] + branch["tm_min"])/2.0
    step_size = (branch["tm_max"] - branch["tm_min"])/(branch["tm_steps"]-1)

    return mid_point + step_size*tm_step
end

function ta_step_to_value(branch, ta_step::Int)
    mid_point = (branch["ta_max"] + branch["ta_min"])/2.0
    step_size = (branch["ta_max"] - branch["ta_min"])/(branch["ta_steps"]-1)

    return mid_point + step_size*ta_step
end

