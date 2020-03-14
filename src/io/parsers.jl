##### Generic Helper Functions #####

function _remove_comment(string)
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
function parse_goc_files(ini_file; scenario_id="")
    files, scenario_id = find_goc_files(ini_file, scenario_id=scenario_id)
    return parse_goc_files(files["con"], files["inl"], files["raw"], files["rop"], ini_file=ini_file, scenario_id=scenario_id)
end

function find_goc_files(ini_file; scenario_id="")
    files = Dict(
        "rop" => "x",
        "raw" => "x",
        "con" => "x",
        "inl" => "x"
    )

    if !endswith(ini_file, ".ini")
        warn(LOGGER, "given init file does not end with .ini, $(ini_file)")
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
                warn(LOGGER, "unknown input given in ini file: $(line)")
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
        info(LOGGER, "no scenario specified, selected directory \"$(scenario_id)\"")
    else
        if !(scenario_id in scenario_dirs)
            error(LOGGER, "$(scenario_id) not found in $(scenario_dirs)")
        end
    end

    for (id, path) in files
        if path == "."
            files[id] = ini_dir
        elseif path == "x"
            files[id] = joinpath(ini_dir, scenario_id)
        else
            error(LOGGER, "unknown file path directive $(path) for file $(id)")
        end
    end

    files["raw"] = joinpath(files["raw"], "case.raw")
    files["rop"] = joinpath(files["rop"], "case.rop")
    files["inl"] = joinpath(files["inl"], "case.inl")
    files["con"] = joinpath(files["con"], "case.con")

    return files, scenario_id
end

function parse_goc_files(con_file, inl_file, raw_file, rop_file; ini_file="", scenario_id="none")
    files = Dict(
        "rop" => rop_file,
        "raw" => raw_file,
        "con" => con_file,
        "inl" => inl_file
    )

    info(LOGGER, "Parsing Files")
    info(LOGGER, "  raw: $(files["raw"])")
    info(LOGGER, "  rop: $(files["rop"])")
    info(LOGGER, "  inl: $(files["inl"])")
    info(LOGGER, "  con: $(files["con"])")

    info(LOGGER, "skipping power models data warnings")
    pm_logger_level = getlevel(getlogger(PowerModels))
    setlevel!(getlogger(PowerModels), "error")
    network_model = PowerModels.parse_file(files["raw"], import_all=true)
    #network_model = parse_psse(files["raw"], import_all=true)
    #@time network_model = parse_psse(files["raw"], import_all=true)

    setlevel!(getlogger(PowerModels), pm_logger_level)

    gen_cost = parse_rop_file(files["rop"])
    response = parse_inl_file(files["inl"])
    contingencies = parse_con_file(files["con"])

    return (ini_file=ini_file, scenario=scenario_id, network=network_model, cost=gen_cost, response=response, contingencies=contingencies, files=files)
end

function parse_goc_opf_files(ini_file; scenario_id="")
    files = Dict(
        "rop" => "x",
        "raw" => "x",
    )

    if !endswith(ini_file, ".ini")
        warn(LOGGER, "given init file does not end with .ini, $(ini_file)")
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
                warn(LOGGER, "unknown input given in ini file: $(line)")
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
        info(LOGGER, "no scenario specified, selected directory \"$(scenario_id)\"")
    else
        if !(scenario_id in scenario_dirs)
            error(LOGGER, "$(scenario_id) not found in $(scenario_dirs)")
        end
    end

    for (id, path) in files
        if path == "."
            files[id] = ini_dir
        elseif path == "x"
            files[id] = joinpath(ini_dir, scenario_id)
        else
            error(LOGGER, "unknown file path directive $(path) for file $(id)")
        end
    end

    files["raw"] = joinpath(files["raw"], "case.raw")
    files["rop"] = joinpath(files["rop"], "case.rop")

    info(LOGGER, "Parsing Files")
    info(LOGGER, "  raw: $(files["raw"])")
    info(LOGGER, "  rop: $(files["rop"])")

    network_model = PowerModels.parse_file(files["raw"], import_all=true)
    gen_cost = parse_rop_file(files["rop"])

    return (ini_file=ini_file, scenario=scenario_id, network=network_model, cost=gen_cost, files=files)
end


##### Unit Inertia and Governor Response Data File Parser (.inl) #####

function parse_inl_file(file::String)
    open(file) do io
        return parse_inl_file(io)
    end
end

function parse_inl_file(io::IO)
    inl_list = []
    for line in readlines(io)
        #line = _remove_comment(line)

        if startswith(strip(line), "0")
            debug(LOGGER, "inl file sentinel found")
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

_rop_sections = [
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

function parse_rop_file(file::String)
    open(file) do io
        return parse_rop_file(io)
    end
end

function parse_rop_file(io::IO)
    active_section_idx = 1
    active_section = _rop_sections[active_section_idx]

    section_data = Dict()
    section_data[active_section.first] = []

    line_idx = 1
    lines = readlines(io)
    while line_idx < length(lines)
        #line = _remove_comment(lines[line_idx])
        line = lines[line_idx]
        if startswith(strip(line), "0")
            debug(LOGGER, "finished reading rop section $(active_section.second) with $(length(section_data[active_section.first])) items")
            active_section_idx += 1
            if active_section_idx > length(_rop_sections)
                debug(LOGGER, "finished reading known rop sections")
                break
            end
            active_section = _rop_sections[active_section_idx]
            section_data[active_section.first] = []
            line_idx += 1
            continue
        end

        if active_section.first == "gen"
            push!(section_data[active_section.first], _parse_rop_gen(line))
        elseif active_section.first == "disptbl"
            push!(section_data[active_section.first], _parse_rop_pg(line))
        elseif active_section.first == "ctbl"
            pwl_line_parts = split(line, ",")
            @assert length(pwl_line_parts) >= 3

            num_pwl_lines = parse(Int, pwl_line_parts[3])
            @assert num_pwl_lines > 0

            pwl_point_lines = lines[line_idx+1:line_idx+num_pwl_lines]
            #pwl_point_lines = remove_comment.(pwl_point_lines)
            push!(section_data[active_section.first], _parse_rop_pwl(pwl_line_parts, pwl_point_lines))
            line_idx += num_pwl_lines
        else
            info(LOGGER, "skipping data line: $(line)")
        end
        line_idx += 1
    end
    return section_data
end

function _parse_rop_gen(line)
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

function _parse_rop_pg(line)
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

function _parse_rop_pwl(pwl_parts, point_lines)
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





##### Contingency Description Data File (.con) #####

# OPEN BRANCH FROM BUS *I TO BUS *J CIRCUIT *1CKT
_branch_contigency_structure = [
    1 => "OPEN",
    2 => "BRANCH",
    3 => "FROM",
    4 => "BUS",
    #5 => "I",
    6 => "TO",
    7 => "BUS",
    #8 => "J",
    9 => "CIRCUIT",
    #10 => "CKT"
]

# REMOVE UNIT *ID FROM BUS *I
_generator_contigency_structure = [
    1 => "REMOVE",
    2 => "UNIT",
    #3 => "ID",
    4 => "FROM",
    5 => "BUS",
    #6 => "I"
]


function parse_con_file(file::String)
    open(file) do io
        return parse_con_file(io)
    end
end

function parse_con_file(io::IO)
    con_lists = []

    tokens = []

    for line in readlines(io)
        #line_tokens = split(strip(_remove_comment(line)))
        line_tokens = split(strip(line))
        #println(line_tokens)
        append!(tokens, line_tokens)
    end

    #println(tokens)

    token_idx = 1
    while token_idx <= length(tokens)
        token = tokens[token_idx]
        if token == "END"
            debug(LOGGER, "end of contingency file found")
            break
        elseif token == "CONTINGENCY"
            # start reading contingencies

            contingency_name = tokens[token_idx+1]
            debug(LOGGER, "reading contingency $(contingency_name)")

            token_idx += 2
            token = tokens[token_idx]
            remaining_tokens = length(tokens) - token_idx

            if token == "OPEN" # branch contingency case
                # OPEN BRANCH FROM BUS *I TO BUS *J CIRCUIT *1CKT

                @assert remaining_tokens >= 9
                branch_tokens = tokens[token_idx:token_idx+9]
                #println(branch_tokens)

                #if !all(branch_tokens[idx] == val for (idx, val) in _branch_contigency_structure) && !all(branch_tokens[idx] == val for (idx, val) in _branch_contigency_structure_alt)
                if any(branch_tokens[idx] != val for (idx, val) in _branch_contigency_structure)
                    error(LOGGER, "incorrect branch contingency structure: $(branch_tokens)")
                end

                bus_i = parse(Int, branch_tokens[5])
                @assert bus_i >= 0

                bus_j = parse(Int, branch_tokens[8])
                @assert bus_j >= 0

                ckt = branch_tokens[10]

                branch_contingency = Dict(
                    "label" => contingency_name,
                    "component" => "branch",
                    "action" => "open",
                    "i" => bus_i,
                    "j" => bus_j,
                    "ckt" => ckt,
                )

                push!(con_lists, branch_contingency)

                token_idx += 9
            elseif token == "REMOVE"
                # REMOVE UNIT *ID FROM BUS *I

                @assert remaining_tokens >= 5
                generator_tokens = tokens[token_idx:token_idx+5]
                #println(generator_tokens)

                if any(generator_tokens[idx] != val for (idx, val) in _generator_contigency_structure)
                    error(LOGGER, "incorrect generator contingency structure: $(generator_tokens)")
                end

                gen_id = generator_tokens[3]

                bus_i = parse(Int, generator_tokens[6])
                @assert bus_i >= 0

                generator_contingency = Dict(
                    "label" => contingency_name,
                    "component" => "generator",
                    "action" => "remove",
                    "id" => gen_id,
                    "i" => bus_i,
                )

                push!(con_lists, generator_contingency)

                token_idx += 5
            elseif token == "END"
                warn(LOGGER, "no action provided for contingency $(contingency_name)")
                token_idx -= 1
            else
                warn(LOGGER, "unrecognized token $(token)")
            end

            token_idx += 1
            token = tokens[token_idx]
            if token != "END"
                error(LOGGER, "expected END token at end of CONTINGENCY, got $(token)")
            end
        else
            warn(LOGGER, "unrecognized token $(token)")
        end
        token_idx += 1
    end

    return con_lists
end




function parse_solution1_file(file::String)
    open(file) do io
        return parse_solution1_file(io)
    end
end

function parse_solution1_file(io::IO)
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
            warn(LOGGER, "skipping blank line in solution1 file ($(idx))")
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
                warn(LOGGER, "skipping line in solution1 file ($(idx)): $(line)")
            end
        end
        idx += 1
    end

    return (bus=bus_data_list, gen=gen_data_list)
end
