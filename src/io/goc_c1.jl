##### GOC C1 Data Tools #####


const C1_VM_EQ_TOL = 1e-4


##### Generic Helper Functions #####

function _remove_psse_comment(string)
    return split(string, "/")[1]
end



##### GOC Initialization File Parser (.ini) #####

"""
Given an `.ini` file and a scenario id (i.e. a directory name), parses the
files from ARPA-e GOC Challenge 1 data format.
* `.raw` network model
* `.rop` generator costs
* `.con` contingency set
* `.inl` generator contingency response parameters
"""
function parse_c1_case(ini_file; scenario_id="")
    files, scenario_id = find_c1_files(ini_file, scenario_id=scenario_id)
    return parse_c1_files(files["con"], files["inl"], files["raw"], files["rop"], ini_file=ini_file, scenario_id=scenario_id)
end

function find_c1_files(ini_file; scenario_id="")
    files = Dict(
        "rop" => "x",
        "raw" => "x",
        "con" => "x",
        "inl" => "x"
    )

    if !endswith(ini_file, ".ini")
        warn(_LOGGER, "given init file does not end with .ini, $(ini_file)")
    end

    open(ini_file) do io
        for line in readlines(io)
            line = strip(line)
            #println(line)
            if startswith(line, "[INPUTS]")
                # do nothing
            elseif startswith(line, "ROP")
                files["rop"] = strip(split(line,"=")[2])
            elseif startswith(line, "RAW")
                files["raw"] = strip(split(line,"=")[2])
            elseif startswith(line, "CON")
                files["con"] = strip(split(line,"=")[2])
            elseif startswith(line, "INL")
                files["inl"] = strip(split(line,"=")[2])
            else
                warn(_LOGGER, "unknown input given in ini file: $(line)")
            end
        end
    end

    #println(files)

    ini_dir = dirname(ini_file)

    #println(ini_dir)
    scenario_dirs = [file for file in readdir(ini_dir) if isdir(joinpath(ini_dir, file))]
    scenario_dirs = sort(scenario_dirs)
    #println(scenario_dirs)

    if length(scenario_id) == 0
        scenario_id = scenario_dirs[1]
        info(_LOGGER, "no scenario specified, selected directory \"$(scenario_id)\"")
    else
        if !(scenario_id in scenario_dirs)
            error(_LOGGER, "$(scenario_id) not found in $(scenario_dirs)")
        end
    end

    for (id, path) in files
        if path == "."
            files[id] = ini_dir
        elseif path == "x"
            files[id] = joinpath(ini_dir, scenario_id)
        else
            error(_LOGGER, "unknown file path directive $(path) for file $(id)")
        end
    end

    files["raw"] = joinpath(files["raw"], "case.raw")
    files["rop"] = joinpath(files["rop"], "case.rop")
    files["inl"] = joinpath(files["inl"], "case.inl")
    files["con"] = joinpath(files["con"], "case.con")

    return files, scenario_id
end

function parse_c1_files(con_file, inl_file, raw_file, rop_file; ini_file="", scenario_id="none")
    files = Dict(
        "rop" => rop_file,
        "raw" => raw_file,
        "con" => con_file,
        "inl" => inl_file
    )

    info(_LOGGER, "Parsing Files")
    info(_LOGGER, "  raw: $(files["raw"])")
    info(_LOGGER, "  rop: $(files["rop"])")
    info(_LOGGER, "  inl: $(files["inl"])")
    info(_LOGGER, "  con: $(files["con"])")

    info(_LOGGER, "skipping power models data warnings")
    pm_logger_level = getlevel(getlogger(PowerModels))
    setlevel!(getlogger(PowerModels), "error")
    network_model = _PM.parse_file(files["raw"], import_all=true)
    #network_model = parse_psse(files["raw"], import_all=true)
    #@time network_model = parse_psse(files["raw"], import_all=true)

    setlevel!(getlogger(PowerModels), pm_logger_level)

    gen_cost = parse_c1_rop_file(files["rop"])
    response = parse_c1_inl_file(files["inl"])
    contingencies = parse_con_file(files["con"])

    return (ini_file=ini_file, scenario=scenario_id, network=network_model, cost=gen_cost, response=response, contingencies=contingencies, files=files)
end

function parse_c1_opf_files(ini_file; scenario_id="")
    files = Dict(
        "rop" => "x",
        "raw" => "x",
    )

    if !endswith(ini_file, ".ini")
        warn(_LOGGER, "given init file does not end with .ini, $(ini_file)")
    end

    open(ini_file) do io
        for line in readlines(io)
            line = strip(line)
            #println(line)
            if startswith(line, "[INPUTS]")
                # do nothing
            elseif startswith(line, "ROP")
                files["rop"] = strip(split(line,"=")[2])
            elseif startswith(line, "RAW")
                files["raw"] = strip(split(line,"=")[2])
            else
                warn(_LOGGER, "unknown input given in ini file: $(line)")
            end
        end
    end

    #println(files)

    ini_dir = dirname(ini_file)

    #println(ini_dir)
    scenario_dirs = [file for file in readdir(ini_dir) if isdir(joinpath(ini_dir, file))]
    scenario_dirs = sort(scenario_dirs)
    #println(scenario_dirs)

    if length(scenario_id) == 0
        scenario_id = scenario_dirs[1]
        info(_LOGGER, "no scenario specified, selected directory \"$(scenario_id)\"")
    else
        if !(scenario_id in scenario_dirs)
            error(_LOGGER, "$(scenario_id) not found in $(scenario_dirs)")
        end
    end

    for (id, path) in files
        if path == "."
            files[id] = ini_dir
        elseif path == "x"
            files[id] = joinpath(ini_dir, scenario_id)
        else
            error(_LOGGER, "unknown file path directive $(path) for file $(id)")
        end
    end

    files["raw"] = joinpath(files["raw"], "case.raw")
    files["rop"] = joinpath(files["rop"], "case.rop")

    info(_LOGGER, "Parsing Files")
    info(_LOGGER, "  raw: $(files["raw"])")
    info(_LOGGER, "  rop: $(files["rop"])")

    network_model = _PM.parse_file(files["raw"], import_all=true)
    gen_cost = parse_c1_rop_file(files["rop"])

    return (ini_file=ini_file, scenario=scenario_id, network=network_model, cost=gen_cost, files=files)
end


##### Unit Inertia and Governor Response Data File Parser (.inl) #####

function parse_c1_inl_file(file::String)
    open(file) do io
        return parse_c1_inl_file(io)
    end
end

function parse_c1_inl_file(io::IO)
    inl_list = []
    for line in readlines(io)
        #line = _remove_psse_comment(line)

        if startswith(strip(line), "0")
            debug(_LOGGER, "inl file sentinel found")
            break
        end
        line_parts = split(line, ",")
        @assert length(line_parts) >= 7

        inl_data = Dict(
            "i"    => parse(Int, line_parts[1]),
            "id"   => strip(line_parts[2]),
            "h"    => strip(line_parts[3]),
            "pmax" => strip(line_parts[4]),
            "pmin" => strip(line_parts[5]),
            "r"    => parse(Float64, line_parts[6]),
            "d"    => strip(line_parts[7])
        )

        @assert inl_data["r"] >= 0.0

        #println(inl_data)
        push!(inl_list, inl_data)
    end
    return inl_list
end




##### Generator Cost Data File Parser (.rop) #####

const _c1_rop_sections = [
    "mod" => "Modification Code",
    "bus_vm" => "Bus Voltage Attributes",
    "shunt_adj" => "Adjustable Bus Shunts",
    "load" => "Bus Loads",
    "load_adj" => "Adjustable Bus Load Tables",
    "gen" => "Generator Dispatch Units",
    "disptbl" => "Active Power Dispatch Tables",
    "gen_reserve" => "Generator Reserve Units",
    "qg" => "Generation Reactive Capability",
    "branch_x" => "Adjustable Branch Reactance",
    "ctbl" => "Piecewise Linear Cost Curve Tables",
    "pwc" => "Piecewise Quadratic Cost Curve Tables",
    "pec" => "Polynomial & Exponential Cost Curve Tables",
    "reserve" => "Period Reserves",
    "branch_flow" => "Branch Flows",
    "int_flow" => "Interface Flows",
    "lin_const" => "Linear Constraint Dependencies",
    "dc_const" => "Two Terminal DC Line Constraint Dependencies",
]

function parse_c1_rop_file(file::String)
    open(file) do io
        return parse_c1_rop_file(io)
    end
end

function parse_c1_rop_file(io::IO)
    active_section_idx = 1
    active_section = _c1_rop_sections[active_section_idx]

    section_data = Dict()
    section_data[active_section.first] = []

    line_idx = 1
    lines = readlines(io)
    while line_idx < length(lines)
        #line = _remove_psse_comment(lines[line_idx])
        line = lines[line_idx]
        if startswith(strip(line), "0")
            debug(_LOGGER, "finished reading rop section $(active_section.second) with $(length(section_data[active_section.first])) items")
            active_section_idx += 1
            if active_section_idx > length(_c1_rop_sections)
                debug(_LOGGER, "finished reading known rop sections")
                break
            end
            active_section = _c1_rop_sections[active_section_idx]
            section_data[active_section.first] = []
            line_idx += 1
            continue
        end

        if active_section.first == "gen"
            push!(section_data[active_section.first], _parse_c1_rop_gen(line))
        elseif active_section.first == "disptbl"
            push!(section_data[active_section.first], _parse_c1_rop_pg(line))
        elseif active_section.first == "ctbl"
            pwl_line_parts = split(line, ",")
            @assert length(pwl_line_parts) >= 3

            num_pwl_lines = parse(Int, pwl_line_parts[3])
            @assert num_pwl_lines > 0

            pwl_point_lines = lines[line_idx+1:line_idx+num_pwl_lines]
            #pwl_point_lines = remove_comment.(pwl_point_lines)
            push!(section_data[active_section.first], _parse_c1_rop_pwl(pwl_line_parts, pwl_point_lines))
            line_idx += num_pwl_lines
        else
            info(_LOGGER, "skipping data line: $(line)")
        end
        line_idx += 1
    end
    return section_data
end

function _parse_c1_rop_gen(line)
    line_parts = split(line, ",")
    @assert length(line_parts) >= 4

    data = Dict(
        "bus"     => parse(Int, line_parts[1]),
        "genid"   => strip(line_parts[2]),
        "disp"    => strip(line_parts[3]),
        "disptbl" => parse(Int, line_parts[4]),
    )

    @assert data["disptbl"] >= 0

    return data
end

function _parse_c1_rop_pg(line)
    line_parts = split(line, ",")
    @assert length(line_parts) >= 7

    data = Dict(
        "tbl"      => parse(Int, line_parts[1]),
        "pmax"     => strip(line_parts[2]),
        "pmin"     => strip(line_parts[3]),
        "fuelcost" => strip(line_parts[4]),
        "ctyp"     => strip(line_parts[5]),
        "status"   => strip(line_parts[6]),
        "ctbl"     => parse(Int, line_parts[7]),
    )

    @assert data["tbl"] >= 0
    @assert data["ctbl"] >= 0

    return data
end

function _parse_c1_rop_pwl(pwl_parts, point_lines)
    @assert length(pwl_parts) >= 2

    points = []

    for point_line in point_lines
        point_line_parts = split(point_line, ",")
        @assert length(point_line_parts) >= 2
        x = parse(Float64, point_line_parts[1])
        y = parse(Float64, point_line_parts[2])

        push!(points, (x=x, y=y))
    end

    data = Dict(
        "ltbl"   =>  parse(Int, pwl_parts[1]),
        "label"  => strip(pwl_parts[2]),
        "points" => points
    )

    @assert data["ltbl"] >= 0

    return data
end






function parse_c1_solution1_file(file::String)
    open(file) do io
        return parse_c1_solution1_file(io)
    end
end

function parse_c1_solution1_file(io::IO)
    bus_data_list = []
    gen_data_list = []

    lines = readlines(io)

    # skip bus list header section
    idx = 1

    separator_count = 0
    skip_next = false

    while idx <= length(lines)
        line = lines[idx]
        if length(strip(line)) == 0
            warn(_LOGGER, "skipping blank line in solution1 file ($(idx))")
        elseif skip_next
            skip_next = false
        elseif startswith(strip(line), "--")
            separator_count += 1
            skip_next = true
        else
            if separator_count == 1
                parts = split(line, ",")
                @assert length(parts) >= 4
                bus_data = (
                    bus = parse(Int, parts[1]),
                    vm = parse(Float64, parts[2]),
                    va = parse(Float64, parts[3]),
                    bcs = parse(Float64, parts[4])
                )
                push!(bus_data_list, bus_data)
            elseif separator_count == 2
                parts = split(line, ",")
                @assert length(parts) >= 4
                gen_data = (
                    bus = parse(Int, parts[1]),
                    id = strip(strip(parts[2]), ['\'', ' ']),
                    pg = parse(Float64, parts[3]),
                    qg = parse(Float64, parts[4])
                )
                push!(gen_data_list, gen_data)
            else
                warn(_LOGGER, "skipping line in solution1 file ($(idx)): $(line)")
            end
        end
        idx += 1
    end

    return (bus=bus_data_list, gen=gen_data_list)
end



##### Transform GOC Data in PM Data #####

"""
transforms files from ARPA-e GOC Challenge 1 data format in to the PowerModels
data format.  This consists of taking the data from multiple data structures
and putting into a network data dictionary.
"""
function build_c1_pm_model(c1_data)
    scenario = c1_data.scenario
    network = c1_data.network

    ##### General Helpers #####
    gen_lookup = Dict(tuple(gen["source_id"][2], strip(gen["source_id"][3])) => gen for (i,gen) in network["gen"])

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



    ##### Link Generator Cost Data #####

    @assert network["per_unit"]
    mva_base = network["baseMVA"]

    dispatch_tbl_lookup = Dict()
    for dispatch_tbl in c1_data.cost["disptbl"]
        dispatch_tbl_lookup[dispatch_tbl["ctbl"]] = dispatch_tbl
    end

    cost_tbl_lookup = Dict()
    for cost_tbl in c1_data.cost["ctbl"]
        cost_tbl_lookup[cost_tbl["ltbl"]] = cost_tbl
    end

    gen_cost_models = Dict()
    for gen_dispatch in c1_data.cost["gen"]
        gen_id = (gen_dispatch["bus"], strip(replace(gen_dispatch["genid"], "'" => "")))
        dispatch_tbl = dispatch_tbl_lookup[gen_dispatch["disptbl"]]
        cost_tbl = cost_tbl_lookup[dispatch_tbl["ctbl"]]

        gen_cost_models[gen_id] = cost_tbl
    end

    if length(gen_cost_models) != length(network["gen"])
        error(_LOGGER, "cost model data missing, network has $(length(network["gen"])) generators, the cost model has $(length(gen_cost_models)) generators")
    end

    for (gen_id, cost_model) in gen_cost_models
        pm_gen = gen_lookup[gen_id]
        pm_gen["model"] = 1
        pm_gen["model_label"] = cost_model["label"]
        pm_gen["ncost"] = length(cost_model["points"])

        #println(cost_model["points"])
        point_list = Float64[]
        for point in cost_model["points"]
            push!(point_list, point.x/mva_base)
            push!(point_list, point.y)
        end
        pm_gen["cost"] = point_list
    end



    ##### Link Generator Participation Data #####

    if length(c1_data.response) != length(network["gen"])
        error(_LOGGER, "generator response model data missing, network has $(length(network["gen"])) generators, the response model has $(length(c1_data.response)) generators")
    end

    for gen_response in c1_data.response
        gen_id = (gen_response["i"], strip(gen_response["id"]))

        pm_gen = gen_lookup[gen_id]

        pm_gen["alpha"] = gen_response["r"]
    end


    ##### Setup Generator Area Group Data #####

    area_gens = Dict{Int,Set{Int}}()
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0 
            gen_bus = network["bus"]["$(gen["gen_bus"])"]
            area = gen_bus["area"]
            if !haskey(area_gens, area)
                area_gens[area] = Set{Int}()
            end
            push!(area_gens[area], gen["index"])
        end
    end
    network["area_gens"] = area_gens


    ##### Flexible Shunt Data #####

    for (i,shunt) in network["shunt"]
        #println(shunt["source_id"])
        if shunt["source_id"][1] == "switched shunt"
            #@assert shunt["source_id"][3] == 0
            @assert shunt["gs"] == 0.0
            shunt["dispatchable"] = true

            bmin = 0.0
            bmax = 0.0
            for (n_name,b_name) in [("n1","b1"),("n2","b2"),("n3","b3"),("n4","b4"),("n5","b5"),("n6","b6"),("n7","b7"),("n8","b8")]
                if shunt[b_name] <= 0.0
                    bmin += shunt[n_name]*shunt[b_name]
                else
                    bmax += shunt[n_name]*shunt[b_name]
                end
            end
            shunt["bmin"] = bmin/mva_base
            shunt["bmax"] = bmax/mva_base
        else
            shunt["dispatchable"] = false
        end
    end



    ##### Add Contingency Lists #####

    generator_ids = []
    branch_ids = []

    for (i,cont) in enumerate(c1_data.contingencies)
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

    # FYI, this breaks output API
    #_PM.propagate_topology_status!(network)

    for (i,shunt) in network["shunt"]
        # test checks if a "switched shunt" in the orginal data model
        if shunt["dispatchable"]
            if shunt["bs"] < shunt["bmin"]
                warn(_LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmin"])")
                shunt["bs"] = shunt["bmin"]
            end
            if shunt["bs"] > shunt["bmax"]
                warn(_LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmax"])")
                shunt["bs"] = shunt["bmax"]
            end
        end
    end

    return network
end


"a simpler version of `build_pm_model` that does not require contingency information"
function build_c1_pm_opf_model(c1_data)
    scenario = c1_data.scenario
    network = c1_data.network

    ##### General Helpers #####

    gen_lookup = Dict(tuple(gen["source_id"][2], strip(gen["source_id"][3])) => gen for (i,gen) in network["gen"])

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



    ##### Link Generator Cost Data #####

    @assert network["per_unit"]
    mva_base = network["baseMVA"]

    dispatch_tbl_lookup = Dict()
    for dispatch_tbl in c1_data.cost["disptbl"]
        dispatch_tbl_lookup[dispatch_tbl["ctbl"]] = dispatch_tbl
    end

    cost_tbl_lookup = Dict()
    for cost_tbl in c1_data.cost["ctbl"]
        cost_tbl_lookup[cost_tbl["ltbl"]] = cost_tbl
    end

    gen_cost_models = Dict()
    for gen_dispatch in c1_data.cost["gen"]
        gen_id = (gen_dispatch["bus"], strip(gen_dispatch["genid"]))
        dispatch_tbl = dispatch_tbl_lookup[gen_dispatch["disptbl"]]
        cost_tbl = cost_tbl_lookup[dispatch_tbl["ctbl"]]

        gen_cost_models[gen_id] = cost_tbl
    end

    if length(gen_cost_models) != length(network["gen"])
        error(_LOGGER, "cost model data missing, network has $(length(network["gen"])) generators, the cost model has $(length(gen_cost_models)) generators")
    end

    for (gen_id, cost_model) in gen_cost_models
        pm_gen = gen_lookup[gen_id]
        pm_gen["model"] = 1
        pm_gen["model_label"] = cost_model["label"]
        pm_gen["ncost"] = length(cost_model["points"])

        #println(cost_model["points"])
        point_list = Float64[]
        for point in cost_model["points"]
            push!(point_list, point.x/mva_base)
            push!(point_list, point.y)
        end
        pm_gen["cost"] = point_list
    end


    ##### Flexible Shunt Data #####

    for (i,shunt) in network["shunt"]
        if shunt["source_id"][1] == "switched shunt"
            @assert shunt["source_id"][3] == 0
            @assert shunt["gs"] == 0.0
            shunt["dispatchable"] = true

            bmin = 0.0
            bmax = 0.0
            for (n_name,b_name) in [("n1","b1"),("n2","b2"),("n3","b3"),("n4","b4"),("n5","b5"),("n6","b6"),("n7","b7"),("n8","b8")]
                if shunt[b_name] <= 0.0
                    bmin += shunt[n_name]*shunt[b_name]
                else
                    bmax += shunt[n_name]*shunt[b_name]
                end
            end
            shunt["bmin"] = bmin/mva_base
            shunt["bmax"] = bmax/mva_base
        else
            shunt["dispatchable"] = false
        end
    end


    ##### Fix Broken Data #####

    _PM.correct_cost_functions!(network)

    # FYI, this breaks output API
    #_PM.propagate_topology_status!(network)

    for (i,shunt) in network["shunt"]
        # test checks if a "switched shunt" in the orginal data model
        if shunt["dispatchable"]
            if shunt["bs"] < shunt["bmin"]
                warn(_LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmin"])")
                shunt["bs"] = shunt["bmin"]
            end
            if shunt["bs"] > shunt["bmax"]
                warn(_LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmax"])")
                shunt["bs"] = shunt["bmax"]
            end
        end
    end


    return network
end




##### Working with GOC Solution Files #####


# this Does not appear to be relivent (03/14/2020)
# function build_contingency_solutions(network, solution)
#     contingency_solutions = Dict{String,Any}()

#     solution = deepcopy(solution)
#     solution["delta"] = 0.0

#     for cont in network["gen_contingencies"]
#         contingency_solutions[cont.label] = solution
#     end

#     for cont in network["branch_contingencies"]
#         contingency_solutions[cont.label] = solution
#     end

#     return contingency_solutions
# end


"checks feasibility criteria of network solution, produces an error if a problem is found"
function check_c1_network_solution(network)
    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            if bus["vm"] > bus["vmax"] || bus["vm"] < bus["vmin"]
                error(_LOGGER, "vm on $(bus["source_id"]) is not in bounds $(bus["vmin"]) to $(bus["vmax"]), given $(bus["vm"])")
            end
        end
    end

    for (i,shunt) in network["shunt"]
        if shunt["status"] != 0
            if haskey(shunt, "dispatchable")
                if shunt["dispatchable"]
                    @assert shunt["gs"] == 0.0
                    @assert haskey(shunt, "bmin") && haskey(shunt, "bmax")
                    if shunt["bs"] > shunt["bmax"] || shunt["bs"] < shunt["bmin"]
                        error(_LOGGER, "bs on $(shunt["source_id"]) is not in bounds $(shunt["bmin"]) to $(shunt["bmax"]), given $(shunt["bs"])")
                    end
                end
            end
        end
    end

    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            if gen["pg"] > gen["pmax"] || gen["pg"] < gen["pmin"]
                error(_LOGGER, "pg on gen $(gen["source_id"]) not in bounds $(gen["pmin"]) to $(gen["pmax"]), given $(gen["pg"])")
            end

            if gen["qg"] > gen["qmax"] || gen["qg"] < gen["qmin"]
                error(_LOGGER, "pg on gen $(gen["source_id"]) not in bounds $(gen["qmin"]) to $(gen["qmax"]), given $(gen["qg"])")
            end
        end
    end
end


"checks feasibility criteria of network solution, corrects when possible"
function correct_c1_solution!(network)

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

    bs_changes = [0.0]
    for (i,shunt) in network["shunt"]
        if haskey(shunt, "dispatchable")
            if shunt["dispatchable"]
                @assert shunt["gs"] == 0.0
                @assert haskey(shunt, "bmin") && haskey(shunt, "bmax")
                if shunt["bs"] > shunt["bmax"]
                    warn(_LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmax"])")
                    push!(bs_changes, shunt["bs"] - shunt["bmax"])
                    shunt["bs"] = shunt["bmax"]
                end
                if shunt["bs"] < shunt["bmin"]
                    warn(_LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmin"])")
                    push!(bs_changes, shunt["bmin"] - shunt["bs"])
                    shunt["bs"] = shunt["bmin"]
                end
            end
        else
            warn(_LOGGER, "shunt $(i) missing dispatchable parameter")
        end
    end

    pg_changes = [0.0]
    qg_changes = [0.0]
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
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

    _c1_summary_changes(network, "base_case", vm_changes, bs_changes, pg_changes, qg_changes)
end


"checks feasibility criteria of contingencies, corrects when possible"
function correct_c1_contingency_solutions!(network, contingency_solutions)
    bus_gens = gens_by_bus(network)

    cont_changes = Int64[]
    cont_vm_changes_max = [0.0]
    cont_bs_changes_max = [0.0]
    cont_pg_changes_max = [0.0]
    cont_qg_changes_max = [0.0]

    for cont_sol in contingency_solutions
        changes = correct_c1_contingency_solution!(network, cont_sol; bus_gens=bus_gens)

        push!(cont_changes, changes.changed)
        push!(cont_vm_changes_max, changes.vm_changes_max)
        push!(cont_bs_changes_max, changes.bs_changes_max)
        push!(cont_pg_changes_max, changes.pg_changes_max)
        push!(cont_qg_changes_max, changes.qg_changes_max)
    end

    println("")

    data = [
        "----",
        "bus",
        "branch",
        "gen_cont",
        "branch_cont",
        "changes_count",
        "vm_max_max",
        "bs_max_max",
        "pg_max_max",
        "qg_max_max",
        "vm_max_mean",
        "bs_max_mean",
        "pg_max_mean",
        "qg_max_mean",
    ]
    println(join(data, ", "))

    data = [
        "DATA_CCS",
        length(network["bus"]),
        length(network["branch"]),
        length(network["gen_contingencies"]),
        length(network["branch_contingencies"]),
        sum(cont_changes),
        maximum(cont_vm_changes_max),
        maximum(cont_bs_changes_max),
        maximum(cont_pg_changes_max),
        maximum(cont_qg_changes_max),
        mean(cont_vm_changes_max),
        mean(cont_bs_changes_max),
        mean(cont_pg_changes_max),
        mean(cont_qg_changes_max),
    ]
    println(join(data, ", "))
end


"checks feasibility criteria of contingencies, corrects when possible"
function correct_c1_contingency_solution!(network, cont_sol; bus_gens = gens_by_bus(network))
    label = cont_sol["label"]
    vm_changes = [0.0]
    for (i,bus) in cont_sol["bus"]
        nw_bus = network["bus"][i]

        if nw_bus["bus_type"] != 4
            if length(bus_gens[i]) > 0
                qg = sum(cont_sol["gen"]["$(gen["index"])"]["qg"] for gen in bus_gens[i])
                qmin = sum(gen["qmin"] for gen in bus_gens[i])
                qmax = sum(gen["qmax"] for gen in bus_gens[i])

                if !isapprox(abs(qmin - qmax), 0.0)
                    if qg >= qmax && bus["vm"] - C1_VM_EQ_TOL/10 > nw_bus["vm"]
                        #println(qmin, " ", qg, " ", qmax)
                        #warn(_LOGGER, "update qg $(qmin) $(qg) $(qmax)")
                        warn(_LOGGER, "update vm on bus $(i) in contingency $(label) to match set-point $(bus["vm"]) -> $(nw_bus["vm"]) due to qg upper bound and vm direction")
                        push!(vm_changes, abs(bus["vm"] - nw_bus["vm"]))
                        bus["vm"] = nw_bus["vm"]
                    end

                    if qg <= qmin && bus["vm"] + C1_VM_EQ_TOL/10 < nw_bus["vm"]
                        #println(qmin, " ", qg, " ", qmax)
                        #warn(_LOGGER, "update qg $(qmin) $(qg) $(qmax)")
                        warn(_LOGGER, "update vm on bus $(i) in contingency $(label) to match set-point $(bus["vm"]) -> $(nw_bus["vm"]) due to qg lower bound and vm direction")
                        push!(vm_changes, abs(bus["vm"] - nw_bus["vm"]))
                        bus["vm"] = nw_bus["vm"]
                    end
                end
            end

            if bus["vm"] > nw_bus["vmax"]
                warn(_LOGGER, "update vm on bus $(i) in contingency $(label) to match ub $(bus["vm"]) -> $(nw_bus["vmax"]) due to out of bounds")
                push!(vm_changes, bus["vm"] - nw_bus["vmax"])
                bus["vm"] = nw_bus["vmax"]
            end

            if bus["vm"] < nw_bus["vmin"]
                warn(_LOGGER, "update vm on bus $(i) in contingency $(label) to match lb $(bus["vm"]) -> $(nw_bus["vmin"]) due to out of bounds")
                push!(vm_changes, nw_bus["vmin"] - bus["vm"])
                bus["vm"] = nw_bus["vmin"]
            end
        else
            bus["vm"] = 0.0
            bus["va"] = 0.0
        end
    end


    bs_changes = [0.0]
    if haskey(cont_sol, "shunt")
        for (i,shunt) in cont_sol["shunt"]
            nw_shunt = network["shunt"][i]
            if haskey(nw_shunt, "dispatchable") && nw_shunt["dispatchable"]
                @assert nw_shunt["gs"] == 0.0
                @assert haskey(nw_shunt, "bmin") && haskey(nw_shunt, "bmax")
                if shunt["bs"] > nw_shunt["bmax"]
                    warn(_LOGGER, "update bs on shunt $(i) in contingency $(label) to be in bounds $(shunt["bs"]) -> $(nw_shunt["bmax"])")
                    push!(bs_changes, shunt["bs"] - nw_shunt["bmax"])
                    shunt["bs"] = nw_shunt["bmax"]
                end
                if shunt["bs"] < nw_shunt["bmin"]
                    warn(_LOGGER, "update bs on shunt $(i) in contingency $(label) to be in bounds $(shunt["bs"]) -> $(nw_shunt["bmin"])")
                    push!(bs_changes, nw_shunt["bmin"] - shunt["bs"])
                    shunt["bs"] = nw_shunt["bmin"]
                end
            end
        end
    end

    response_gens = Set{Int}()

    gen_id = -1
    if cont_sol["cont_type"] == "gen"
        gen_id = cont_sol["cont_comp_id"]
        nw_gen = network["gen"]["$(gen_id)"]
        nw_gen_bus = network["bus"]["$(nw_gen["gen_bus"])"]
        response_gens = network["area_gens"][nw_gen_bus["area"]]
    end

    if cont_sol["cont_type"] == "branch"
        branch_id = cont_sol["cont_comp_id"]
        nw_branch = network["branch"]["$(branch_id)"]
        nw_fr_bus = network["bus"]["$(nw_branch["f_bus"])"]
        nw_to_bus = network["bus"]["$(nw_branch["t_bus"])"]
        response_gens = Set()
        if haskey(network["area_gens"], nw_fr_bus["area"])
            response_gens = network["area_gens"][nw_fr_bus["area"]]
        end
        if nw_fr_bus["area"] != nw_to_bus["area"] && haskey(network["area_gens"], nw_to_bus["area"])
            response_gens = union(response_gens, network["area_gens"][nw_to_bus["area"]])
        end
    end

    pg_changes = [0.0]
    qg_changes = [0.0]
    delta = cont_sol["delta"]
    for (i,gen) in cont_sol["gen"]
        nw_gen = network["gen"][i]

        if !(nw_gen["gen_status"] == 0 || (gen_id >= 0 && nw_gen["index"] == gen_id))
            bus_id = nw_gen["gen_bus"]
            nw_bus = network["bus"]["$(bus_id)"]

            if gen["qg"] < nw_gen["qmax"] && gen["qg"] > nw_gen["qmin"]
                bus = cont_sol["bus"]["$(bus_id)"]
                #if !isapprox(bus["vm"], nw_bus["vm"])
                if !isapprox(bus["vm"], nw_bus["vm"], atol=C1_VM_EQ_TOL/2)
                    #debug(_LOGGER, "bus $(bus_id) : vm_base $(nw_bus["vm"]) - vm $(bus["vm"]) : reactive bounds $(nw_gen["qmin"]) - $(gen["qg"]) - $(nw_gen["qmax"])")
                    warn(_LOGGER, "update vm on bus $(bus_id) in contingency $(label) to match base case $(bus["vm"]) -> $(nw_bus["vm"]) due to within reactive bounds")
                end
                bus["vm"] = nw_bus["vm"]
            end

            pg_calc = nw_gen["pg"] 
            if nw_gen["index"] in response_gens
                pg_calc += nw_gen["alpha"]*delta
                pg_calc = max(pg_calc, nw_gen["pmin"])
                pg_calc = min(pg_calc, nw_gen["pmax"])
            end

            if !isapprox(gen["pg"], pg_calc, atol=1e-5)
                warn(_LOGGER, "pg value on gen $(i) $(nw_gen["source_id"]) in contingency $(label) is not consistent with the computed value given:$(gen["pg"]) calc:$(pg_calc)")
            end

            if gen["pg"] > nw_gen["pmax"]
                warn(_LOGGER, "update pg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match ub $(gen["pg"]) -> $(nw_gen["pmax"])")
                push!(pg_changes, gen["pg"] - nw_gen["pmax"])
                gen["pg"] = nw_gen["pmax"]
            end

            if gen["pg"] < nw_gen["pmin"]
                warn(_LOGGER, "update pg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match lb $(gen["pg"]) -> $(nw_gen["pmin"])")
                push!(pg_changes, nw_gen["pmin"] - gen["pg"])
                gen["pg"] = nw_gen["pmin"]
            end

            if gen["qg"] > nw_gen["qmax"]
                warn(_LOGGER, "update qg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match ub $(gen["qg"]) -> $(nw_gen["qmax"])")
                push!(qg_changes, gen["qg"] - nw_gen["qmax"])
                gen["qg"] = nw_gen["qmax"]
            end

            if gen["qg"] < nw_gen["qmin"]
                warn(_LOGGER, "update qg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match lb $(gen["qg"]) -> $(nw_gen["qmin"])")
                push!(qg_changes, nw_gen["qmin"] - gen["qg"])
                gen["qg"] = nw_gen["qmin"]
            end
        else
            gen["pg"] = 0.0
            gen["qg"] = 0.0
        end
    end

    # test imbalance on branch conts after flow corrections
    #=
    network_tmp = deepcopy(network)
    cont_type = cont_sol["cont_type"]
    cont_idx = cont_sol["cont_comp_id"]
    network_tmp["branch"]["$(cont_idx)"]["br_status"] = 0
    _PM.update_data!(network_tmp, cont_sol)
    deltas = calc_c1_power_balance_deltas!(network_tmp)
    println("$(label): $(deltas)")
    =#

    cont_changed = length(vm_changes) > 1 || length(bs_changes) > 1 || length(pg_changes) > 1 || length(qg_changes) > 1

    if cont_changed
        _c1_summary_changes(network, label, vm_changes, bs_changes, pg_changes, qg_changes)
    end

    return (changed=Int(cont_changed), vm_changes_max=maximum(vm_changes), bs_changes_max=maximum(bs_changes), pg_changes_max=maximum(pg_changes), qg_changes_max=maximum(qg_changes))
end


function _c1_summary_changes(network, contingency, vm_changes, bs_changes, pg_changes, qg_changes)
    #println("")

    data = [
        "----",
        "contingency",
        "bus",
        "branch",
        "gen_cont",
        "branch_cont",
        "vm_count",
        "bs_count",
        "pg_count",
        "qg_count",
        "vm_max",
        "bs_max",
        "pg_max",
        "qg_max",
        "vm_mean",
        "bs_mean",
        "pg_mean",
        "qg_mean",
    ]
    info(_LOGGER, join(data, ", "))

    data = [
        "DATA_CHANGES",
        contingency,
        length(network["bus"]),
        length(network["branch"]),
        length(network["gen_contingencies"]),
        length(network["branch_contingencies"]),
        length(vm_changes)-1,
        length(bs_changes)-1,
        length(pg_changes)-1,
        length(qg_changes)-1,
        maximum(vm_changes),
        maximum(bs_changes),
        maximum(pg_changes),
        maximum(qg_changes),
        mean(vm_changes),
        mean(bs_changes),
        mean(pg_changes),
        mean(qg_changes),
    ]
    info(_LOGGER, join(data, ", "))
end



function write_c1_solution(c1_data, network, contingencies; solution_1="", solution_2="")
    error(_LOGGER, "the write_solution function has been replaced by write_solution1, write_solution2")
end


function write_c1_solution1(network; output_dir="", solution_file="solution1.txt")
    if length(output_dir) > 0
        solution_path = joinpath(output_dir, solution_file)
    else
        solution_path = solution_file
    end

    open(solution_path, "w") do sol1
        base_mva = network["baseMVA"]

        bus_switched_shunt_b = Dict(i => 0.0 for (i,bus) in network["bus"])
        for (i,shunt) in network["shunt"]
            # test checks if a "switched shunt" in the orginal data model
            if shunt["dispatchable"] && shunt["status"] == 1
                @assert shunt["gs"] == 0.0
                bus_switched_shunt_b["$(shunt["shunt_bus"])"] += shunt["bs"]
            end
        end

        write(sol1, "-- bus section\n")
        write(sol1, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
        for (i,bus) in network["bus"]
            write(sol1, "$(bus["index"]), $(bus["vm"]), $(rad2deg(bus["va"])), $(base_mva*bus_switched_shunt_b[i])\n")
        end

        write(sol1, "-- generator section\n")
        write(sol1, "i, id, p(MW), q(MVAR)\n")
        for (i,gen) in network["gen"]
            bus_index = gen["source_id"][2]
            gen_id = gen["source_id"][3]
            write(sol1, "$(bus_index), $(gen_id), $(base_mva*gen["pg"]), $(base_mva*gen["qg"])\n")
        end
    end
end


function write_c1_solution2(network, contingencies; output_dir="", solution_file="solution2.txt")
    if length(output_dir) > 0
        solution_path = joinpath(output_dir, solution_file)
    else
        solution_path = solution_file
    end

    open(solution_path, "w") do sol2
        base_mva = network["baseMVA"]

        for cont_solution in contingencies
            write_c1_solution2_contingency(sol2, network, cont_solution)
        end
    end

    return solution_path
end


function write_c1_solution2_contingency(io::IO, network, contingency_solution)
    base_mva = network["baseMVA"]

    bus_switched_shunt_b = Dict(i => 0.0 for (i,bus) in network["bus"])
    if haskey(network, "shunt") && haskey(contingency_solution, "shunt")
        for (i,nw_shunt) in network["shunt"]
            if nw_shunt["dispatchable"] && nw_shunt["status"] == 1
                #@assert nw_shunt["gs"] == 0.0
                shunt = contingency_solution["shunt"][i]
                bus_switched_shunt_b["$(nw_shunt["shunt_bus"])"] += shunt["bs"]
            end
        end
    end

    write(io, "-- contingency\n")
    write(io, "label\n")
    write(io, "$(contingency_solution["label"])\n")

    write(io, "-- bus section\n")
    write(io, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
    for (i,nw_bus) in network["bus"]
        vm = 0.0
        va = 0.0

        if haskey(contingency_solution["bus"], i)
            bus = contingency_solution["bus"][i]
            vm = bus["vm"]
            va = bus["va"]
        end

        write(io, "$(nw_bus["index"]), $(vm), $(rad2deg(va)), $(base_mva*bus_switched_shunt_b[i])\n")
    end

    write(io, "-- generator section\n")
    write(io, "i, id, p(MW), q(MVAR)\n")
    for (i,nw_gen) in network["gen"]

        nw_gen = network["gen"][i]
        bus_index = nw_gen["source_id"][2]
        gen_id = nw_gen["source_id"][3]

        pg = 0.0
        qg = 0.0

        if haskey(contingency_solution["gen"], i)
            gen = contingency_solution["gen"][i]
            pg = gen["pg"]
            qg = gen["qg"]
        end

        write(io, "$(bus_index), $(gen_id), $(base_mva*pg), $(base_mva*qg)\n")
    end

    write(io, "-- delta section\n")
    write(io, "delta(MW)\n")
    write(io, "$(base_mva*contingency_solution["delta"])\n")
end


function c1_combine_files(files, output_file_name; output_dir="")
    if length(output_dir) > 0
        output_path = joinpath(output_dir, output_file_name)
    else
        output_path = output_file_name
    end

    open(output_path, "w") do output
        for file in files
            open(file, "r") do input
                for line in readlines(input, keep=true)
                    write(output, line)
                end
            end
        end
    end

    return output_path
end


function remove_c1_files(files)
    for file in files
        if isfile(file)
            info(_LOGGER, "deleting: $(file)")
            rm(file)
        else
            info(_LOGGER, "skipping file: $(file)")
        end
    end
end


function remove_c1_solution_files(;output_dir="", solution1_file="solution1.txt", solution2_file="solution2.txt")
    if length(output_dir) > 0
        solution1_path = joinpath(output_dir, solution1_file)
    else
        solution1_path = solution1_file
    end

    if isfile(solution1_path)
        info(_LOGGER, "deleting: $(solution1_path)")
        rm(solution1_path)
    else
        info(_LOGGER, "skipping file: $(solution1_path)")
    end


    if length(output_dir) > 0
        solution2_path = joinpath(output_dir, solution2_file)
    else
        solution2_path = solution2_file
    end

    if isfile(solution2_path)
        info(_LOGGER, "deleting: $(solution2_path)")
        rm(solution2_path)
    else
        info(_LOGGER, "skipping file: $(solution2_path)")
    end
end


function remove_c1_detail_file(;output_dir="", detail_file="detail.csv")
    if length(output_dir) > 0
        detail_file_path = joinpath(output_dir, detail_file)
    else
        detail_file_path = detail_file
    end

    if isfile(detail_file_path)
        info(_LOGGER, "deleting: $(detail_file_path)")
        rm(detail_file_path)
    else
        info(_LOGGER, "skipping file: $(detail_file_path)")
    end
end


function write_c1_file_paths(files; output_dir="", solution1_file="solution1.txt", solution2_file="solution2.txt")
    if length(output_dir) > 0
        solution1_path = joinpath(output_dir, solution1_file)
    else
        solution1_path = solution1_file
    end

    if length(output_dir) > 0
        solution2_path = joinpath(output_dir, solution2_file)
    else
        solution2_path = solution2_file
    end

    files["sol1"] = solution1_path
    files["sol2"] = solution2_path
    open("files.json", "w") do input
        JSON.print(input, files)
    end
end



function read_c1_solution1(network; output_dir="", state_file="solution1.txt")
    if length(output_dir) > 0
        solution1_path = joinpath(output_dir, state_file)
    else
        solution1_path = state_file
    end

    return build_c1_pm_solution(network, solution1_path)
end

function build_c1_pm_solution(network, goc_sol_file::String)
    info(_LOGGER, "loading solution file: $(goc_sol_file)")
    goc_sol = parse_c1_solution1_file(goc_sol_file)

    info(_LOGGER, "converting GOC solution to PowerModels solution")
    pm_sol = build_c1_pm_solution(network, goc_sol)

    return pm_sol
end

function build_c1_pm_solution(network, goc_sol)
    bus_lookup = Dict(parse(Int, bus["source_id"][2]) => bus for (i,bus) in network["bus"])
    gen_lookup = Dict((gen["source_id"][2], strip(gen["source_id"][3])) => gen for (i,gen) in network["gen"])
    shunt_lookup = Dict{Int,Any}()
    for (i,shunt) in network["shunt"]
        if shunt["source_id"][1] == "switched shunt"
            @assert shunt["source_id"][3] == 0
            shunt_lookup[shunt["source_id"][2]] = shunt
        end
    end

    base_mva = network["baseMVA"]

    bus_data = Dict{String,Any}()
    shunt_data = Dict{String,Any}()
    for bus_sol in goc_sol.bus
        pm_bus = bus_lookup[bus_sol.bus]
        bus_data["$(pm_bus["index"])"] = Dict(
            "vm" => bus_sol.vm,
            "va" => deg2rad(bus_sol.va)
        )

        if haskey(shunt_lookup, bus_sol.bus)
            pm_shunt = shunt_lookup[bus_sol.bus]
            shunt_data["$(pm_shunt["index"])"] = Dict(
                "gs" => 0.0,
                "bs" => bus_sol.bcs/base_mva
            )
        else
            @assert bus_sol.bcs == 0.0
        end
    end

    gen_data = Dict{String,Any}()
    for gen_sol in goc_sol.gen
        pm_gen = gen_lookup[(gen_sol.bus, gen_sol.id)]
        gen_data["$(pm_gen["index"])"] = Dict(
            "pg" => gen_sol.pg/base_mva,
            "qg" => gen_sol.qg/base_mva
        )
    end

    solution = Dict(
        "per_unit" => true,
        "bus" => bus_data,
        "shunt" => shunt_data,
        "gen" => gen_data
    )

    return solution
end


