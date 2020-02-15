##### Generic Helper Functions #####

function remove_comment(string)
    return split(string, "/")[1]
end



##### GOC Initialization File Parser (.ini) #####

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
        #line = remove_comment(line)

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

rop_sections = [
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
    active_section = rop_sections[active_section_idx]

    section_data = Dict()
    section_data[active_section.first] = []

    line_idx = 1
    lines = readlines(io)
    while line_idx < length(lines)
        #line = remove_comment(lines[line_idx])
        line = lines[line_idx]
        if startswith(strip(line), "0")
            debug(LOGGER, "finished reading rop section $(active_section.second) with $(length(section_data[active_section.first])) items")
            active_section_idx += 1
            if active_section_idx > length(rop_sections)
                debug(LOGGER, "finished reading known rop sections")
                break
            end
            active_section = rop_sections[active_section_idx]
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
branch_contigency_structure = [
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
generator_contigency_structure = [
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
        #line_tokens = split(strip(remove_comment(line)))
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

                #if !all(branch_tokens[idx] == val for (idx, val) in branch_contigency_structure) && !all(branch_tokens[idx] == val for (idx, val) in branch_contigency_structure_alt)
                if any(branch_tokens[idx] != val for (idx, val) in branch_contigency_structure)
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

                if any(generator_tokens[idx] != val for (idx, val) in generator_contigency_structure)
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








# ****** PSSE Parser is below this line ******


_LOGGER = LOGGER


"""
    _init_bus!(bus, id)

Initializes a `bus` of id `id` with default values given in the PSS(R)E
specification.
"""
function _init_bus!(bus::Dict{String,Any}, id::Int)
    bus["bus_i"] = id
    bus["bus_type"] = 1
    bus["area"] = 1
    bus["vm"] = 1.0
    bus["va"] = 0.0
    bus["base_kv"] = 1.0
    bus["zone"] = 1
    bus["name"] = "            "
    bus["vmax"] = 1.1
    bus["vmin"] = 0.9
    bus["index"] = id
end


"""
    _get_bus_value(bus_i, field, pm_data)

Returns the value of `field` of `bus_i` from the PowerModels data. Requires
"bus" Dict to already be populated.
"""
function _get_bus_value(bus_i, field, pm_data)
    if isa(pm_data["bus"], Array)
        for bus in pm_data["bus"]
            if bus["index"] == bus_i
                return bus[field]
            end
        end
    elseif isa(pm_data["bus"], Dict)
        for (k, bus) in pm_data["bus"]
            if bus["index"] == bus_i
                return bus[field]
            end
        end
    end

    Memento.warn(_LOGGER, "Could not find bus $bus_i, returning 0 for field $field")
    return 0
end


"""
    _find_max_bus_id(pm_data)

Returns the maximum bus id in `pm_data`
"""
function _find_max_bus_id(pm_data::Dict)::Int
    max_id = 0
    for bus in pm_data["bus"]
        if bus["index"] > max_id && !endswith(bus["name"], "starbus")
            max_id = bus["index"]
        end
    end

    return max_id
end


"""
    create_starbus(pm_data, transformer)

Creates a starbus from a given three-winding `transformer`. "source_id" is given
by `["bus_i", "name", "I", "J", "K", "CKT"]` where "bus_i" and "name" are the
modified names for the starbus, and "I", "J", "K" and "CKT" come from the
originating transformer, in the PSS(R)E transformer specification.
"""
function _create_starbus_from_transformer(pm_data::Dict, transformer::Dict)::Dict
    starbus = Dict{String,Any}()

    # transformer starbus ids will be one order of magnitude larger than highest real bus id
    base = convert(Int, 10 ^ ceil(log10(abs(_find_max_bus_id(pm_data)))))
    starbus_id = transformer["I"] + base

    _init_bus!(starbus, starbus_id)

    starbus["name"] = "$(transformer["I"]) starbus"

    starbus["vm"] = transformer["VMSTAR"]
    starbus["va"] = transformer["ANSTAR"]
    starbus["bus_type"] = transformer["STAT"]
    starbus["area"] = _get_bus_value(transformer["I"], "area", pm_data)
    starbus["zone"] = _get_bus_value(transformer["I"], "zone", pm_data)
    starbus["source_id"] = push!(["transformer", starbus["bus_i"], starbus["name"]], transformer["I"], transformer["J"], transformer["K"], transformer["CKT"])

    return starbus
end


"Imports remaining keys from `data_in` into `data_out`, excluding keys in `exclude`"
function _import_remaining!(data_out::Dict, data_in::Dict, import_all::Bool; exclude=[])
    if import_all
        for (k, v) in data_in
            if !(k in exclude)
                if isa(v, Array)
                    for (n, item) in enumerate(v)
                        if isa(item, Dict)
                            _import_remaining!(item, item, import_all)
                            if !("index" in keys(item))
                                item["index"] = n
                            end
                        end
                    end
                elseif isa(v, Dict)
                    _import_remaining!(v, v, import_all)
                end
                data_out[lowercase(k)] = v
                delete!(data_in, k)
            end
        end
    end
end


"""
    _psse2pm_branch!(pm_data, pti_data)

Parses PSS(R)E-style Branch data into a PowerModels-style Dict. "source_id" is
given by `["I", "J", "CKT"]` in PSS(R)E Branch specification.
"""
function _psse2pm_branch!(pm_data::Dict, pti_data::Dict, import_all::Bool)


    pm_data["branch"] = []
    if haskey(pti_data, "BRANCH")
        for (i, branch) in enumerate(pti_data["BRANCH"])
            sub_data = Dict{String,Any}()

            sub_data["f_bus"] = pop!(branch, "I")
            sub_data["t_bus"] = pop!(branch, "J")
            sub_data["br_r"] = pop!(branch, "R")
            sub_data["br_x"] = pop!(branch, "X")
            sub_data["g_fr"] = pop!(branch, "GI")
            sub_data["b_fr"] = branch["BI"] == 0. && branch["B"] != 0. ? branch["B"] / 2 : pop!(branch, "BI")
            sub_data["g_to"] = pop!(branch, "GJ")
            sub_data["b_to"] = branch["BJ"] == 0. && branch["B"] != 0. ? branch["B"] / 2 : pop!(branch, "BJ")
            sub_data["rate_a"] = pop!(branch, "RATEA")
            sub_data["rate_b"] = pop!(branch, "RATEB")
            sub_data["rate_c"] = pop!(branch, "RATEC")
            sub_data["tap"] = 1.0
            sub_data["shift"] = 0.0
            sub_data["br_status"] = pop!(branch, "ST")
            sub_data["angmin"] = 0.0
            sub_data["angmax"] = 0.0
            sub_data["transformer"] = false

            sub_data["source_id"] = ["branch", sub_data["f_bus"], sub_data["t_bus"], pop!(branch, "CKT")]
            sub_data["index"] = i

            _import_remaining!(sub_data, branch, import_all; exclude=["B", "BI", "BJ"])

            if sub_data["rate_a"] == 0.0
                delete!(sub_data, "rate_a")
            end
            if sub_data["rate_b"] == 0.0
                delete!(sub_data, "rate_b")
            end
            if sub_data["rate_c"] == 0.0
                delete!(sub_data, "rate_c")
            end

            push!(pm_data["branch"], sub_data)
        end
    end
end


"""
    _psse2pm_generator!(pm_data, pti_data)

Parses PSS(R)E-style Generator data in a PowerModels-style Dict. "source_id" is
given by `["I", "ID"]` in PSS(R)E Generator specification.
"""
function _psse2pm_generator!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    pm_data["gen"] = []
    if haskey(pti_data, "GENERATOR")
        for gen in pti_data["GENERATOR"]
            sub_data = Dict{String,Any}()

            sub_data["gen_bus"] = pop!(gen, "I")
            sub_data["gen_status"] = pop!(gen, "STAT")
            sub_data["pg"] = pop!(gen, "PG")
            sub_data["qg"] = pop!(gen, "QG")
            sub_data["vg"] = pop!(gen, "VS")
            sub_data["mbase"] = pop!(gen, "MBASE")
            sub_data["pmin"] = pop!(gen, "PB")
            sub_data["pmax"] = pop!(gen, "PT")
            sub_data["qmin"] = pop!(gen, "QB")
            sub_data["qmax"] = pop!(gen, "QT")

            # Default Cost functions
            sub_data["model"] = 2
            sub_data["startup"] = 0.0
            sub_data["shutdown"] = 0.0
            sub_data["ncost"] = 2
            sub_data["cost"] = [1.0, 0.0]

            sub_data["source_id"] = ["generator", sub_data["gen_bus"], pop!(gen, "ID")]
            sub_data["index"] = length(pm_data["gen"]) + 1

            _import_remaining!(sub_data, gen, import_all)

            push!(pm_data["gen"], sub_data)
        end
    end
end


"""
    _psse2pm_bus!(pm_data, pti_data)

Parses PSS(R)E-style Bus data into a PowerModels-style Dict. "source_id" is given
by ["I", "NAME"] in PSS(R)E Bus specification.
"""
function _psse2pm_bus!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    pm_data["bus"] = []
    if haskey(pti_data, "BUS")
        for bus in pti_data["BUS"]
            sub_data = Dict{String,Any}()

            sub_data["bus_i"] = bus["I"]
            sub_data["bus_type"] = pop!(bus, "IDE")
            sub_data["area"] = pop!(bus, "AREA")
            sub_data["vm"] = pop!(bus, "VM")
            sub_data["va"] = pop!(bus, "VA")
            sub_data["base_kv"] = pop!(bus, "BASKV")
            sub_data["zone"] = pop!(bus, "ZONE")
            sub_data["name"] = pop!(bus, "NAME")
            sub_data["vmax"] = pop!(bus, "NVHI")
            sub_data["vmin"] = pop!(bus, "NVLO")

            sub_data["source_id"] = ["bus", "$(bus["I"])"]
            sub_data["index"] = pop!(bus, "I")

            _import_remaining!(sub_data, bus, import_all)

            push!(pm_data["bus"], sub_data)
        end
    end
end


"""
    _psse2pm_load!(pm_data, pti_data)

Parses PSS(R)E-style Load data into a PowerModels-style Dict. "source_id" is given
by `["I", "ID"]` in the PSS(R)E Load specification.
"""
function _psse2pm_load!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    pm_data["load"] = []
    if haskey(pti_data, "LOAD")
        for load in pti_data["LOAD"]
            sub_data = Dict{String,Any}()

            sub_data["load_bus"] = pop!(load, "I")
            sub_data["pd"] = pop!(load, "PL")
            sub_data["qd"] = pop!(load, "QL")
            sub_data["status"] = pop!(load, "STATUS")

            sub_data["source_id"] = ["load", sub_data["load_bus"], pop!(load, "ID")]
            sub_data["index"] = length(pm_data["load"]) + 1

            _import_remaining!(sub_data, load, import_all)

            push!(pm_data["load"], sub_data)
        end
    end
end


"""
    _psse2pm_shunt!(pm_data, pti_data)

Parses PSS(R)E-style Fixed and Switched Shunt data into a PowerModels-style
Dict. "source_id" is given by `["I", "ID"]` for Fixed Shunts, and `["I", "SWREM"]`
for Switched Shunts, as given by the PSS(R)E Fixed and Switched Shunts
specifications.
"""
function _psse2pm_shunt!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    pm_data["shunt"] = []

    if haskey(pti_data, "FIXED SHUNT")
        for shunt in pti_data["FIXED SHUNT"]
            sub_data = Dict{String,Any}()

            sub_data["shunt_bus"] = pop!(shunt, "I")
            sub_data["gs"] = pop!(shunt, "GL")
            sub_data["bs"] = pop!(shunt, "BL")
            sub_data["status"] = pop!(shunt, "STATUS")

            sub_data["source_id"] = ["fixed shunt", sub_data["shunt_bus"], pop!(shunt, "ID")]
            sub_data["index"] = length(pm_data["shunt"]) + 1

            _import_remaining!(sub_data, shunt, import_all)

            push!(pm_data["shunt"], sub_data)
        end
    end

    if haskey(pti_data, "SWITCHED SHUNT")
        Memento.info(_LOGGER, "Switched shunt converted to fixed shunt, with default value gs=0.0")

        for shunt in pti_data["SWITCHED SHUNT"]
            sub_data = Dict{String,Any}()

            sub_data["shunt_bus"] = pop!(shunt, "I")
            sub_data["gs"] = 0.0
            sub_data["bs"] = pop!(shunt, "BINIT")
            sub_data["status"] = pop!(shunt, "STAT")

            sub_data["source_id"] = ["switched shunt", sub_data["shunt_bus"], pop!(shunt, "SWREM")]
            sub_data["index"] = length(pm_data["shunt"]) + 1

            _import_remaining!(sub_data, shunt, import_all)

            push!(pm_data["shunt"], sub_data)
        end
    end
end


"""
    _psse2pm_transformer!(pm_data, pti_data)

Parses PSS(R)E-style Transformer data into a PowerModels-style Dict. "source_id"
is given by `["I", "J", "K", "CKT", "winding"]`, where "winding" is 0 if
transformer is two-winding, and 1, 2, or 3 for three-winding, and the remaining
keys are defined in the PSS(R)E Transformer specification.
"""
function _psse2pm_transformer!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    if !haskey(pm_data, "branch")
        pm_data["branch"] = []
    end

    if haskey(pti_data, "TRANSFORMER")
        for transformer in pti_data["TRANSFORMER"]
            if transformer["K"] == 0  # Two-winding Transformers
                sub_data = Dict{String,Any}()

                sub_data["f_bus"] = transformer["I"]
                sub_data["t_bus"] = transformer["J"]

                # Unit Transformations
                if transformer["CZ"] == 1  # "for resistance and reactance in pu on system MVA base and winding voltage base"
                    br_r, br_x = transformer["R1-2"], transformer["X1-2"]
                else  # NOT "for resistance and reactance in pu on system MVA base and winding voltage base"
                    if transformer["CZ"] == 3  # "for transformer load loss in watts and impedance magnitude in pu on a specified MVA base and winding voltage base."
                        br_r = 1e-6 * transformer["R1-2"] / transformer["SBASE1-2"]
                        br_x = sqrt(transformer["X1-2"]^2 - br_r^2)
                    else
                        br_r, br_x = transformer["R1-2"], transformer["X1-2"]
                    end
                    br_r *= (transformer["NOMV1"]^2 / _get_bus_value(transformer["I"], "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE1-2"])
                    br_x *= (transformer["NOMV1"]^2 / _get_bus_value(transformer["I"], "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE1-2"])
                end

                sub_data["br_r"] = br_r
                sub_data["br_x"] = br_x

                sub_data["g_fr"] = pop!(transformer, "MAG1")
                sub_data["b_fr"] = pop!(transformer, "MAG2")
                sub_data["g_to"] = 0.0
                sub_data["b_to"] = 0.0

                sub_data["rate_a"] = pop!(transformer, "RATA1")
                sub_data["rate_b"] = pop!(transformer, "RATB1")
                sub_data["rate_c"] = pop!(transformer, "RATC1")

                if sub_data["rate_a"] == 0.0
                    delete!(sub_data, "rate_a")
                end
                if sub_data["rate_b"] == 0.0
                    delete!(sub_data, "rate_b")
                end
                if sub_data["rate_c"] == 0.0
                    delete!(sub_data, "rate_c")
                end

                sub_data["tap"] = pop!(transformer, "WINDV1") / pop!(transformer, "WINDV2")
                sub_data["shift"] = pop!(transformer, "ANG1")

                # Unit Transformations
                if transformer["CW"] != 1  # NOT "for off-nominal turns ratio in pu of winding bus base voltage"
                    sub_data["tap"] *= _get_bus_value(transformer["J"], "base_kv", pm_data) / _get_bus_value(transformer["I"], "base_kv", pm_data)
                    if transformer["CW"] == 3  # "for off-nominal turns ratio in pu of nominal winding voltage, NOMV1, NOMV2 and NOMV3."
                        sub_data["tap"] *= transformer["NOMV1"] / transformer["NOMV2"]
                    end
                end

                sub_data["br_status"] = transformer["STAT"]

                sub_data["angmin"] = 0.0
                sub_data["angmax"] = 0.0

                sub_data["source_id"] = ["transformer", pop!(transformer, "I"), pop!(transformer, "J"), pop!(transformer, "K"), pop!(transformer, "CKT"), 0]
                sub_data["transformer"] = true
                sub_data["index"] = length(pm_data["branch"]) + 1

                _import_remaining!(sub_data, transformer, import_all; exclude=["I", "J", "K", "CZ", "CW", "R1-2", "R2-3", "R3-1",
                                                                              "X1-2", "X2-3", "X3-1", "SBASE1-2", "SBASE2-3",
                                                                              "SBASE3-1", "MAG1", "MAG2", "STAT", "NOMV1", "NOMV2"])

                push!(pm_data["branch"], sub_data)
            else  # Three-winding Transformers
                bus_id1, bus_id2, bus_id3 = transformer["I"], transformer["J"], transformer["K"]

                # Creates a starbus (or "dummy" bus) to which each winding of the transformer will connect
                starbus = _create_starbus_from_transformer(pm_data, transformer)
                push!(pm_data["bus"], starbus)

                # Create 3 branches from a three winding transformer (one for each winding, which will each connect to the starbus)
                br_r12, br_r23, br_r31 = transformer["R1-2"], transformer["R2-3"], transformer["R3-1"]
                br_x12, br_x23, br_x31 = transformer["X1-2"], transformer["X2-3"], transformer["X3-1"]

                # Unit Transformations
                if transformer["CZ"] == 3  # "for transformer load loss in watts and impedance magnitude in pu on a specified MVA base and winding voltage base."
                    br_r12 *= 1e-6 / transformer["SBASE1-2"]
                    br_r23 *= 1e-6 / transformer["SBASE2-3"]
                    br_r31 *= 1e-6 / transformer["SBASE3-1"]

                    br_x12 = sqrt(br_x12^2 - br_r12^2)
                    br_x23 = sqrt(br_x23^2 - br_r23^2)
                    br_x31 = sqrt(br_x31^2 - br_r31^2)
                end

                # Unit Transformations
                if transformer["CZ"] != 1  # NOT "for resistance and reactance in pu on system MVA base and winding voltage base"
                    br_r12 *= (transformer["NOMV1"]^2 / _get_bus_value(bus_id1, "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE1-2"])
                    br_r23 *= (transformer["NOMV2"]^2 / _get_bus_value(bus_id2, "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE2-3"])
                    br_r31 *= (transformer["NOMV3"]^2 / _get_bus_value(bus_id3, "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE3-1"])

                    br_x12 *= (transformer["NOMV1"]^2 / _get_bus_value(bus_id1, "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE1-2"])
                    br_x23 *= (transformer["NOMV2"]^2 / _get_bus_value(bus_id2, "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE2-3"])
                    br_x31 *= (transformer["NOMV3"]^2 / _get_bus_value(bus_id3, "base_kv", pm_data)^2) * (pm_data["baseMVA"] / transformer["SBASE3-1"])
                end

                # See "Power System Stability and Control", ISBN: 0-07-035958-X, Eq. 6.72
                Zr_p = 1/2 * (br_r12 - br_r23 + br_r31)
                Zr_s = 1/2 * (br_r23 - br_r31 + br_r12)
                Zr_t = 1/2 * (br_r31 - br_r12 + br_r23)
                Zx_p = 1/2 * (br_x12 - br_x23 + br_x31)
                Zx_s = 1/2 * (br_x23 - br_x31 + br_x12)
                Zx_t = 1/2 * (br_x31 - br_x12 + br_x23)

                # Build each of the three transformer branches
                for (m, (bus_id, br_r, br_x)) in enumerate(zip([bus_id1, bus_id2, bus_id3], [Zr_p, Zr_s, Zr_t], [Zx_p, Zx_s, Zx_t]))
                    sub_data = Dict{String,Any}()

                    sub_data["f_bus"] = bus_id
                    sub_data["t_bus"] = starbus["bus_i"]

                    sub_data["br_r"] = br_r
                    sub_data["br_x"] = br_x

                    sub_data["g_fr"] = m == 1 ? pop!(transformer, "MAG1") : 0.0
                    sub_data["b_fr"] = m == 1 ? pop!(transformer, "MAG2") : 0.0
                    sub_data["g_to"] = 0.0
                    sub_data["b_to"] = 0.0

                    sub_data["rate_a"] = pop!(transformer, "RATA$m")
                    sub_data["rate_b"] = pop!(transformer, "RATB$m")
                    sub_data["rate_c"] = pop!(transformer, "RATC$m")

                    if sub_data["rate_a"] == 0.0
                        delete!(sub_data, "rate_a")
                    end
                    if sub_data["rate_b"] == 0.0
                        delete!(sub_data, "rate_b")
                    end
                    if sub_data["rate_c"] == 0.0
                        delete!(sub_data, "rate_c")
                    end

                    sub_data["tap"] = pop!(transformer, "WINDV$m")
                    sub_data["shift"] = pop!(transformer, "ANG$m")

                    # Unit Transformations
                    if transformer["CW"] != 1  # NOT "for off-nominal turns ratio in pu of winding bus base voltage"
                        sub_data["tap"] /= _get_bus_value(bus_id, "base_kv", pm_data)
                        if transformer["CW"] == 3  # "for off-nominal turns ratio in pu of nominal winding voltage, NOMV1, NOMV2 and NOMV3."
                            sub_data["tap"] *= transformer["NOMV$m"]
                        end
                    end

                    sub_data["br_status"] = transformer["STAT"]

                    sub_data["angmin"] = 0.0
                    sub_data["angmax"] = 0.0

                    sub_data["source_id"] = ["transformer", transformer["I"], transformer["J"], transformer["K"], transformer["CKT"], m]
                    sub_data["transformer"] = true
                    sub_data["index"] = length(pm_data["branch"]) + 1

                    _import_remaining!(sub_data, transformer, import_all; exclude=["I", "J", "K", "CZ", "CW", "R1-2", "R2-3", "R3-1",
                                                                                  "X1-2", "X2-3", "X3-1", "SBASE1-2", "SBASE2-3", "CKT",
                                                                                  "SBASE3-1", "MAG1", "MAG2", "STAT","NOMV1", "NOMV2",
                                                                                  "NOMV3", "WINDV1", "WINDV2", "WINDV3", "RATA1",
                                                                                  "RATA2", "RATA3", "RATB1", "RATB2", "RATB3", "RATC1",
                                                                                  "RATC2", "RATC3", "ANG1", "ANG2", "ANG3"])

                    push!(pm_data["branch"], sub_data)
                end
            end
        end
    end
end




"""
    _psse2pm_dcline!(pm_data, pti_data)

Parses PSS(R)E-style Two-Terminal and VSC DC Lines data into a PowerModels
compatible Dict structure by first converting them to a simple DC Line Model.
For Two-Terminal DC lines, "source_id" is given by `["IPR", "IPI", "NAME"]` in the
PSS(R)E Two-Terminal DC specification. For Voltage Source Converters, "source_id"
is given by `["IBUS1", "IBUS2", "NAME"]`, where "IBUS1" is "IBUS" of the first
converter bus, and "IBUS2" is the "IBUS" of the second converter bus, in the
PSS(R)E Voltage Source Converter specification.
"""
function _psse2pm_dcline!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    pm_data["dcline"] = []

    if haskey(pti_data, "TWO-TERMINAL DC")
        for dcline in pti_data["TWO-TERMINAL DC"]
            Memento.info(_LOGGER, "Two-Terminal DC lines are supported via a simple *lossless* dc line model approximated by two generators.")
            sub_data = Dict{String,Any}()

            # Unit conversions?
            power_demand = dcline["MDC"] == 1 ? abs(dcline["SETVL"]) : dcline["MDC"] == 2 ? abs(dcline["SETVL"] / pop!(dcline, "VSCHD") / 1000) : 0

            sub_data["f_bus"] = dcline["IPR"]
            sub_data["t_bus"] = dcline["IPI"]
            sub_data["br_status"] = pop!(dcline, "MDC") == 0 ? 0 : 1
            sub_data["pf"] = power_demand
            sub_data["pt"] = power_demand
            sub_data["qf"] = 0.0
            sub_data["qt"] = 0.0
            sub_data["vf"] = _get_bus_value(pop!(dcline, "IPR"), "vm", pm_data)
            sub_data["vt"] = _get_bus_value(pop!(dcline, "IPI"), "vm", pm_data)

            sub_data["pminf"] = 0.0
            sub_data["pmaxf"] = dcline["SETVL"] > 0 ? power_demand : -power_demand
            sub_data["pmint"] = pop!(dcline, "SETVL") > 0 ? -power_demand : power_demand
            sub_data["pmaxt"] = 0.0

            anmn = []
            for key in ["ANMNR", "ANMNI"]
                if abs(dcline[key]) <= 90.
                    push!(anmn, pop!(dcline, key))
                else
                    push!(anmn, 0)
                    Memento.warn(_LOGGER, "$key outside reasonable limits, setting to 0 degress")
                end
            end

            sub_data["qmaxf"] = 0.0
            sub_data["qmaxt"] = 0.0
            sub_data["qminf"] = -max(abs(sub_data["pminf"]), abs(sub_data["pmaxf"])) * cosd(anmn[1])
            sub_data["qmint"] = -max(abs(sub_data["pmint"]), abs(sub_data["pmaxt"])) * cosd(anmn[2])

            # Can we use "number of bridges in series (NBR/NBI)" to compute a loss?
            sub_data["loss0"] = 0.0
            sub_data["loss1"] = 0.0

            # Costs (set to default values)
            sub_data["startup"] = 0.0
            sub_data["shutdown"] = 0.0
            sub_data["ncost"] = 3
            sub_data["cost"] = [0.0, 0.0, 0.0]
            sub_data["model"] = 2

            sub_data["source_id"] = ["two-terminal dc", sub_data["f_bus"], sub_data["t_bus"], pop!(dcline, "NAME")]
            sub_data["index"] = length(pm_data["dcline"]) + 1

            _import_remaining!(sub_data, dcline, import_all)

            push!(pm_data["dcline"], sub_data)
        end
    end

    if haskey(pti_data, "VOLTAGE SOURCE CONVERTER")
        Memento.info(_LOGGER, "VSC-HVDC lines are supported via a dc line model approximated by two generators and an associated loss.")
        for dcline in pti_data["VOLTAGE SOURCE CONVERTER"]
            # Converter buses : is the distinction between ac and dc side meaningful?
            dcside, acside = dcline["CONVERTER BUSES"]

            # PowerWorld conversion from PTI to matpower seems to create two
            # artificial generators from a VSC, but it is not clear to me how
            # the value of "pg" is determined and adds shunt to the DC-side bus.
            sub_data = Dict{String,Any}()

            # VSC intended to be one or bi-directional?
            sub_data["f_bus"] = pop!(dcside, "IBUS")
            sub_data["t_bus"] = pop!(acside, "IBUS")
            sub_data["br_status"] = pop!(dcline, "MDC") == 0 || pop!(dcside, "TYPE") == 0 || pop!(acside, "TYPE") == 0 ? 0 : 1

            sub_data["pf"] = 0.0
            sub_data["pt"] = 0.0

            sub_data["qf"] = 0.0
            sub_data["qt"] = 0.0

            sub_data["vf"] = pop!(dcside, "MODE") == 1 ? pop!(dcside, "ACSET") : 1.0
            sub_data["vt"] = pop!(acside, "MODE") == 1 ? pop!(acside, "ACSET") : 1.0

            sub_data["pmaxf"] = dcside["SMAX"] == 0.0 && dcside["IMAX"] == 0.0 ? max(abs(dcside["MAXQ"]), abs(dcside["MINQ"])) : min(pop!(dcside, "IMAX"), pop!(dcside, "SMAX"))
            sub_data["pmaxt"] = acside["SMAX"] == 0.0 && acside["IMAX"] == 0.0 ? max(abs(acside["MAXQ"]), abs(acside["MINQ"])) : min(pop!(acside, "IMAX"), pop!(acside, "SMAX"))
            sub_data["pminf"] = -sub_data["pmaxf"]
            sub_data["pmint"] = -sub_data["pmaxt"]

            sub_data["qminf"] = pop!(dcside, "MINQ")
            sub_data["qmaxf"] = pop!(dcside, "MAXQ")
            sub_data["qmint"] = pop!(acside, "MINQ")
            sub_data["qmaxt"] = pop!(acside, "MAXQ")

            sub_data["loss0"] = (pop!(dcside, "ALOSS") + pop!(acside, "ALOSS") + pop!(dcside, "MINLOSS") + pop!(acside, "MINLOSS")) * 1e-3
            sub_data["loss1"] = (pop!(dcside, "BLOSS") + pop!(acside, "BLOSS")) * 1e-3 # how to include resistance?

            # Costs (set to default values)
            sub_data["startup"] = 0.0
            sub_data["shutdown"] = 0.0
            sub_data["ncost"] = 3
            sub_data["cost"] = [0.0, 0.0, 0.0]
            sub_data["model"] = 2

            sub_data["source_id"] = ["vsc dc", sub_data["f_bus"], sub_data["t_bus"], pop!(dcline, "NAME")]
            sub_data["index"] = length(pm_data["dcline"]) + 1

            _import_remaining!(sub_data, dcline, import_all)

            push!(pm_data["dcline"], sub_data)
        end
    end
end


function _psse2pm_storage!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    pm_data["storage"] = []
end

function _psse2pm_switch!(pm_data::Dict, pti_data::Dict, import_all::Bool)
    pm_data["switch"] = []
end


"""
    _pti_to_powermodels!(pti_data)

Converts PSS(R)E-style data parsed from a PTI raw file, passed by `pti_data`
into a format suitable for use internally in PowerModels. Imports all remaining
data from the PTI file if `import_all` is true (Default: false).
"""
function _pti_to_powermodels!(pti_data::Dict; import_all=false, validate=true)::Dict
    pm_data = Dict{String,Any}()

    rev = pop!(pti_data["CASE IDENTIFICATION"][1], "REV")

    pm_data["per_unit"] = false
    pm_data["source_type"] = "pti"
    pm_data["source_version"] = "$rev"
    pm_data["baseMVA"] = pop!(pti_data["CASE IDENTIFICATION"][1], "SBASE")
    pm_data["name"] = pop!(pti_data["CASE IDENTIFICATION"][1], "NAME")

    _import_remaining!(pm_data, pti_data["CASE IDENTIFICATION"][1], import_all)

    _psse2pm_bus!(pm_data, pti_data, import_all)
    _psse2pm_load!(pm_data, pti_data, import_all)
    _psse2pm_shunt!(pm_data, pti_data, import_all)
    _psse2pm_generator!(pm_data, pti_data, import_all)
    _psse2pm_branch!(pm_data, pti_data, import_all)
    _psse2pm_transformer!(pm_data, pti_data, import_all)
    _psse2pm_dcline!(pm_data, pti_data, import_all)
    _psse2pm_storage!(pm_data, pti_data, import_all)
    _psse2pm_switch!(pm_data, pti_data, import_all)

    _import_remaining!(pm_data, pti_data, import_all; exclude=[
        "CASE IDENTIFICATION", "BUS", "LOAD", "FIXED SHUNT",
        "SWITCHED SHUNT", "GENERATOR","BRANCH", "TRANSFORMER",
        "TWO-TERMINAL DC", "VOLTAGE SOURCE CONVERTER"
    ])

    # update lookup structure
    for (k, v) in pm_data
        if isa(v, Array)
            #println("updating $(k)")
            dict = Dict{String,Any}()
            for item in v
                @assert("index" in keys(item))
                dict[string(item["index"])] = item
            end
            pm_data[k] = dict
        end
    end

    if validate
        correct_network_data!(pm_data)
    end

    return pm_data
end


"Parses directly from file"
function parse_psse(filename::String; kwargs...)::Dict
    pm_data = open(filename) do f
        parse_psse(f; kwargs...)
    end

    return pm_data
end


"Parses directly from iostream"
function parse_psse(io::IO; kwargs...)::Dict
    pti_data = parse_pti(io)
    #@time pti_data = parse_pti(io)
    pm = _pti_to_powermodels!(pti_data; kwargs...)
    #@time pm = _pti_to_powermodels!(pti_data; kwargs...)
    return pm
end








"""
A list of data file sections in the order that they appear in a PTI v33 file
"""
const _pti_sections = ["CASE IDENTIFICATION", "BUS", "LOAD", "FIXED SHUNT",
    "GENERATOR", "BRANCH", "TRANSFORMER", "AREA INTERCHANGE",
    "TWO-TERMINAL DC", "VOLTAGE SOURCE CONVERTER", "IMPEDANCE CORRECTION",
    "MULTI-TERMINAL DC", "MULTI-SECTION LINE", "ZONE", "INTER-AREA TRANSFER",
    "OWNER", "FACTS CONTROL DEVICE", "SWITCHED SHUNT", "GNE DEVICE",
    "INDUCTION MACHINE"]


const _transaction_dtypes = [("IC", Int64), ("SBASE", Float64), ("REV", Int64),
    ("XFRRAT", Float64), ("NXFRAT", Float64), ("BASFRQ", Float64)]

const _bus_dtypes = [("I", Int64), ("NAME", String), ("BASKV", Float64),
    ("IDE", Int64), ("AREA", Int64), ("ZONE", Int64), ("OWNER", Int64),
    ("VM", Float64), ("VA", Float64), ("NVHI", Float64), ("NVLO", Float64),
    ("EVHI", Float64), ("EVLO", Float64)]

const _load_dtypes = [("I", Int64), ("ID", String), ("STATUS", Int64),
    ("AREA", Int64), ("ZONE", Int64), ("PL", Float64), ("QL", Float64),
    ("IP", Float64), ("IQ", Float64), ("YP", Float64), ("YQ", Float64),
    ("OWNER", Int64), ("SCALE", Int64), ("INTRPT", Int64)]

const _fixed_shunt_dtypes = [("I", Int64), ("ID", String), ("STATUS", Int64),
    ("GL", Float64), ("BL", Float64)]

const _generator_dtypes = [("I", Int64), ("ID", String), ("PG", Float64),
    ("QG", Float64), ("QT", Float64), ("QB", Float64), ("VS", Float64),
    ("IREG", Int64), ("MBASE", Float64), ("ZR", Float64), ("ZX", Float64),
    ("RT", Float64), ("XT", Float64), ("GTAP", Float64), ("STAT", Int64),
    ("RMPCT", Float64), ("PT", Float64), ("PB", Float64), ("O1", Int64),
    ("F1", Float64), ("O2", Int64), ("F2", Float64), ("O3", Int64),
    ("F3", Float64), ("O4", Int64), ("F4", Float64), ("WMOD", Int64),
    ("WPF", Float64)]

const _branch_dtypes = [("I", Int64), ("J", Int64), ("CKT", String),
    ("R", Float64), ("X", Float64), ("B", Float64), ("RATEA", Float64),
    ("RATEB", Float64), ("RATEC", Float64), ("GI", Float64), ("BI", Float64),
    ("GJ", Float64), ("BJ", Float64), ("ST", Int64), ("MET", Int64),
    ("LEN", Float64), ("O1", Int64), ("F1", Float64), ("O2", Int64),
    ("F2", Float64), ("O3", Int64), ("F3", Float64), ("O4", Int64),
    ("F4", Float64)]

const _transformer_dtypes = [("I", Int64), ("J", Int64), ("K", Int64),
    ("CKT", String), ("CW", Int64), ("CZ", Int64), ("CM", Int64),
    ("MAG1", Float64), ("MAG2", Float64), ("NMETR", Int64), ("NAME", String),
    ("STAT", Int64), ("O1", Int64), ("F1", Float64), ("O2", Int64),
    ("F2", Float64), ("O3", Int64), ("F3", Float64), ("O4", Int64),
    ("F4", Float64), ("VECGRP", String)]

const _transformer_3_1_dtypes = [("R1-2", Float64), ("X1-2", Float64),
    ("SBASE1-2", Float64), ("R2-3", Float64), ("X2-3", Float64),
    ("SBASE2-3", Float64), ("R3-1", Float64), ("X3-1", Float64),
    ("SBASE3-1", Float64), ("VMSTAR", Float64), ("ANSTAR", Float64)]

const _transformer_3_2_dtypes = [("WINDV1", Float64), ("NOMV1", Float64),
    ("ANG1", Float64), ("RATA1", Float64), ("RATB1", Float64),
    ("RATC1", Float64), ("COD1", Int64), ("CONT1", Int64), ("RMA1", Float64),
    ("RMI1", Float64), ("VMA1", Float64), ("VMI1", Float64), ("NTP1", Float64),
    ("TAB1", Int64), ("CR1", Float64), ("CX1", Float64), ("CNXA1", Float64)]

const _transformer_3_3_dtypes = [("WINDV2", Float64), ("NOMV2", Float64),
    ("ANG2", Float64), ("RATA2", Float64), ("RATB2", Float64),
    ("RATC2", Float64), ("COD2", Int64), ("CONT2", Int64), ("RMA2", Float64),
    ("RMI2", Float64), ("VMA2", Float64), ("VMI2", Float64), ("NTP2", Float64),
    ("TAB2", Int64), ("CR2", Float64), ("CX2", Float64), ("CNXA2", Float64)]

const _transformer_3_4_dtypes = [("WINDV3", Float64), ("NOMV3", Float64),
    ("ANG3", Float64), ("RATA3", Float64), ("RATB3", Float64),
    ("RATC3", Float64), ("COD3", Int64), ("CONT3", Int64), ("RMA3", Float64),
    ("RMI3", Float64), ("VMA3", Float64), ("VMI3", Float64), ("NTP3", Float64),
    ("TAB3", Int64), ("CR3", Float64), ("CX3", Float64), ("CNXA3", Float64)]

const _transformer_2_1_dtypes = [("R1-2", Float64), ("X1-2", Float64), 
    ("SBASE1-2", Float64)]

const _transformer_2_2_dtypes = [("WINDV1", Float64), ("NOMV1", Float64),
    ("ANG1", Float64), ("RATA1", Float64), ("RATB1", Float64),
    ("RATC1", Float64), ("COD1", Int64), ("CONT1", Int64), ("RMA1", Float64),
    ("RMI1", Float64), ("VMA1", Float64), ("VMI1", Float64), ("NTP1", Float64),
    ("TAB1", Int64), ("CR1", Float64), ("CX1", Float64), ("CNXA1", Float64)]

const _transformer_2_3_dtypes = [("WINDV2", Float64), ("NOMV2", Float64)]

const _area_interchange_dtypes = [("I", Int64), ("ISW", Int64), 
    ("PDES", Float64), ("PTOL", Float64), ("ARNAME", String)]

const _two_terminal_line_dtypes = [("NAME", String), ("MDC", Int64),
    ("RDC", Float64), ("SETVL", Float64), ("VSCHD", Float64),
    ("VCMOD", Float64), ("RCOMP", Float64), ("DELTI", Float64),
    ("METER", String), ("DCVMIN", Float64), ("CCCITMX", Int64),
    ("CCCACC", Float64), ("IPR", Int64), ("NBR", Int64), ("ANMXR", Float64),
    ("ANMNR", Float64), ("RCR", Float64), ("XCR", Float64), ("EBASR", Float64),
    ("TRR", Float64), ("TAPR", Float64), ("TMXR", Float64), ("TMNR", Float64),
    ("STPR", Float64), ("ICR", Int64), ("IFR", Int64), ("ITR", Int64),
    ("IDR", String), ("XCAPR", Float64), ("IPI", Int64), ("NBI", Int64),
    ("ANMXI", Float64), ("ANMNI", Float64), ("RCI", Float64), ("XCI", Float64),
    ("EBASI", Float64), ("TRI", Float64), ("TAPI", Float64), ("TMXI", Float64),
    ("TMNI", Float64), ("STPI", Float64), ("ICI", Int64), ("IFI", Int64),
    ("ITI", Int64), ("IDI", String), ("XCAPI", Float64)]

const _vsc_line_dtypes = [("NAME", String), ("MDC", Int64), ("RDC", Float64),
    ("O1", Int64), ("F1", Float64), ("O2", Int64), ("F2", Float64),
    ("O3", Int64), ("F3", Float64), ("O4", Int64), ("F4", Float64)]

const _vsc_subline_dtypes = [("IBUS", Int64), ("TYPE", Int64), ("MODE", Int64),
    ("DCSET", Float64), ("ACSET", Float64), ("ALOSS", Float64),
    ("BLOSS", Float64), ("MINLOSS", Float64), ("SMAX", Float64),
    ("IMAX", Float64), ("PWF", Float64), ("MAXQ", Float64), ("MINQ", Float64),
    ("REMOT", Int64), ("RMPCT", Float64)]

const _impedance_correction_dtypes = [("I", Int64), ("T1", Float64),
    ("F1", Float64), ("T2", Float64), ("F2", Float64), ("T3", Float64),
    ("F3", Float64), ("T4", Float64), ("F4", Float64), ("T5", Float64),
    ("F5", Float64), ("T6", Float64), ("F6", Float64), ("T7", Float64),
    ("F7", Float64), ("T8", Float64), ("F8", Float64), ("T9", Float64),
    ("F9", Float64), ("T10", Float64), ("F10", Float64), ("T11", Float64),
    ("F11", Float64)]

const _multi_term_main_dtypes = [("NAME", String), ("NCONV", Int64),
    ("NDCBS", Int64), ("NDCLN", Int64), ("MDC", Int64), ("VCONV", Int64),
    ("VCMOD", Float64), ("VCONVN", Float64)]

const _multi_term_nconv_dtypes = [("IB", Int64), ("N", Int64),
    ("ANGMX", Float64), ("ANGMN", Float64), ("RC", Float64), ("XC", Float64),
    ("EBAS", Float64), ("TR", Float64), ("TAP", Float64), ("TPMX", Float64),
    ("TPMN", Float64), ("TSTP", Float64), ("SETVL", Float64),
    ("DCPF", Float64), ("MARG", Float64), ("CNVCOD", Int64)]

const _multi_term_ndcbs_dtypes = [("IDC", Int64), ("IB", Int64),
    ("AREA", Int64), ("ZONE", Int64), ("DCNAME", String), ("IDC2", Int64),
    ("RGRND", Float64), ("OWNER", Int64)]

const _multi_term_ndcln_dtypes = [("IDC", Int64), ("JDC", Int64),
    ("DCCKT", String), ("MET", Int64), ("RDC", Float64), ("LDC", Float64)]

const _multi_section_dtypes = [("I", Int64), ("J", Int64), ("ID", String),
    ("MET", Int64), ("DUM1", Int64), ("DUM2", Int64), ("DUM3", Int64),
    ("DUM4", Int64), ("DUM5", Int64), ("DUM6", Int64), ("DUM7", Int64),
    ("DUM8", Int64), ("DUM9", Int64)]

const _zone_dtypes = [("I", Int64), ("ZONAME", String)]

const _interarea_dtypes = [("ARFROM", Int64), ("ARTO", Int64),
    ("TRID", String), ("PTRAN", Float64)]

const _owner_dtypes = [("I", Int64), ("OWNAME", String)]

const _FACTS_dtypes = [("NAME", String), ("I", Int64), ("J", Int64),
    ("MODE", Int64), ("PDES", Float64), ("QDES", Float64), ("VSET", Float64),
    ("SHMX", Float64), ("TRMX", Float64), ("VTMN", Float64), ("VTMX", Float64),
    ("VSMX", Float64), ("IMX", Float64), ("LINX", Float64), ("RMPCT", Float64),
    ("OWNER", Int64), ("SET1", Float64), ("SET2", Float64), ("VSREF", Int64),
    ("REMOT", Int64), ("MNAME", String)]

const _switched_shunt_dtypes = [("I", Int64), ("MODSW", Int64),
    ("ADJM", Int64), ("STAT", Int64), ("VSWHI", Float64), ("VSWLO", Float64),
    ("SWREM", Int64), ("RMPCT", Float64), ("RMIDNT", String), ("BINIT", Float64),
    ("N1", Int64), ("B1", Float64), ("N2", Int64), ("B2", Float64),
    ("N3", Int64), ("B3", Float64), ("N4", Int64), ("B4", Float64),
    ("N5", Int64), ("B5", Float64), ("N6", Int64), ("B6", Float64),
    ("N7", Int64), ("B7", Float64), ("N8", Int64), ("B8", Float64)]

# TODO: Account for multiple lines in GNE Device entries
const _gne_device_dtypes = [("NAME", String), ("MODEL", String),
    ("NTERM", Int64), ("BUSi", Int64), ("NREAL", Int64), ("NINTG", Int64),
    ("NCHAR", Int64), ("STATUS", Int64), ("OWNER", Int64), ("NMETR", Int64),
    ("REALi", Float64), ("INTGi", Int64), ("CHARi", String)]

const _induction_machine_dtypes = [("I", Int64), ("ID", String),
    ("STAT", Int64), ("SCODE", Int64), ("DCODE", Int64), ("AREA", Int64),
    ("ZONE", Int64), ("OWNER", Int64), ("TCODE", Int64), ("BCODE", Int64),
    ("MBASE", Float64), ("RATEKV", Float64), ("PCODE", Int64),
    ("PSET", Float64), ("H", Float64), ("A", Float64), ("B", Float64),
    ("D", Float64), ("E", Float64), ("RA", Float64), ("XA", Float64),
    ("XM", Float64), ("R1", Float64), ("X1", Float64), ("R2", Float64),
    ("X2", Float64), ("X3", Float64), ("E1", Float64), ("SE1", Float64),
    ("E2", Float64), ("SE2", Float64), ("IA1", Float64), ("IA2", Float64),
    ("XAMULT", Float64)]


"""
lookup array of data types for PTI file sections given by
`field_name`, as enumerated by PSS/E Program Operation Manual.
"""
const _pti_dtypes = Dict{String,Array}(
    "BUS" => _bus_dtypes,
    "LOAD" => _load_dtypes,
    "FIXED SHUNT" => _fixed_shunt_dtypes,
    "GENERATOR" => _generator_dtypes,
    "BRANCH" => _branch_dtypes,
    "TRANSFORMER" => _transformer_dtypes,
    "TRANSFORMER TWO-WINDING LINE 1" => _transformer_2_1_dtypes,
    "TRANSFORMER TWO-WINDING LINE 2" => _transformer_2_2_dtypes,
    "TRANSFORMER TWO-WINDING LINE 3" => _transformer_2_3_dtypes,
    "TRANSFORMER THREE-WINDING LINE 1" => _transformer_3_1_dtypes,
    "TRANSFORMER THREE-WINDING LINE 2" => _transformer_3_2_dtypes,
    "TRANSFORMER THREE-WINDING LINE 3" => _transformer_3_3_dtypes,
    "TRANSFORMER THREE-WINDING LINE 4" => _transformer_3_4_dtypes,
    "AREA INTERCHANGE" => _area_interchange_dtypes,
    "TWO-TERMINAL DC" => _two_terminal_line_dtypes,
    "VOLTAGE SOURCE CONVERTER" => _vsc_line_dtypes,
    "VOLTAGE SOURCE CONVERTER SUBLINES" => _vsc_subline_dtypes,
    "IMPEDANCE CORRECTION" => _impedance_correction_dtypes,
    "MULTI-TERMINAL DC" => _multi_term_main_dtypes,
    "MULTI-TERMINAL DC NCONV" => _multi_term_nconv_dtypes,
    "MULTI-TERMINAL DC NDCBS" => _multi_term_ndcbs_dtypes,
    "MULTI-TERMINAL DC NDCLN" => _multi_term_ndcln_dtypes,
    "MULTI-SECTION LINE" => _multi_section_dtypes,
    "ZONE" => _zone_dtypes,
    "INTER-AREA TRANSFER" => _interarea_dtypes,
    "OWNER" => _owner_dtypes,
    "FACTS CONTROL DEVICE" => _FACTS_dtypes,
    "SWITCHED SHUNT" => _switched_shunt_dtypes,
    "CASE IDENTIFICATION" => _transaction_dtypes,
    "GNE DEVICE" => _gne_device_dtypes,
    "INDUCTION MACHINE" => _induction_machine_dtypes
)



const _default_case_identification = Dict("IC" => 0, "SBASE" => 100.0,
    "REV" => 33, "XFRRAT" => 0, "NXFRAT" => 0, "BASFRQ" => 60)

const _default_bus = Dict("BASKV" => 0.0, "IDE" => 1, "AREA" => 1, "ZONE" => 1,
    "OWNER" => 1, "VM" => 1.0, "VA" => 0.0, "NVHI" => 1.1, "NVLO" => 0.9,
    "EVHI" => 1.1, "EVLO" => 0.9, "NAME" => "            ")

const _default_load = Dict("ID" => 1, "STATUS" => 1, "PL" => 0.0, "QL" => 0.0,
    "IP" => 0.0, "IQ" => 0.0, "YP" => 0.0, "YQ" => 0.0, "SCALE" => 1,
    "INTRPT" => 0, "AREA" => nothing, "ZONE" => nothing, "OWNER" => nothing)

const _default_fixed_shunt = Dict("ID" => 1, "STATUS" => 1, "GL" => 0.0, "BL" => 0.0)

const _default_generator = Dict("ID" => 1, "PG" => 0.0, "QG" => 0.0, "QT" => 9999.0,
    "QB" => -9999.0, "VS" => 1.0, "IREG" => 0, "MBASE" => nothing, "ZR" => 0.0,
    "ZX" => 1.0, "RT" => 0.0, "XT" => 0.0, "GTAP" => 1.0, "STAT" => 1,
    "RMPCT" => 100.0, "PT" => 9999.0, "PB" => -9999.0, "O1" => nothing,
    "O2" => 0, "O3" => 0, "O4" => 0, "F1" => 1.0,"F2" => 1.0, "F3" => 1.0,
    "F4" => 1.0, "WMOD" => 0, "WPF" => 1.0)

const _default_branch = Dict("CKT" => 1, "B" => 0.0, "RATEA" => 0.0,
    "RATEB" => 0.0, "RATEC" => 0.0, "GI" => 0.0, "BI" => 0.0, "GJ" => 0.0,
    "BJ" => 0.0, "ST" => 1, "MET" => 1, "LEN" => 0.0, "O1" => nothing,
    "O2" => 0, "O3" => 0, "O4" => 0, "F1" => 1.0, "F2" => 1.0, "F3" => 1.0,
    "F4" => 1.0)

const _default_transformer = Dict("K" => 0, "CKT" => 1, "CW" => 1, "CZ" => 1,
    "CM" => 1, "MAG1" => 0.0, "MAG2" => 0.0, "NMETR" => 2,
    "NAME" => "            ", "STAT" => 1, "O1" => nothing, "O2" => 0,
    "O3" => 0, "O4" => 0, "F1" => 1.0, "F2" => 1.0, "F3" => 1.0, "F4" => 1.0,
    "VECGRP" => "            ", "R1-2" => 0.0, "SBASE1-2" => nothing,
    "R2-3" => 0.0, "SBASE2-3" => nothing, "R3-1" => 0.0, "SBASE3-1" => nothing,
    "VMSTAR" => 1.0, "ANSTAR" => 0.0,
    "WINDV1" => nothing,
    "NOMV1" => 0.0, "ANG1" => 0.0, "RATA1" => 0.0,
    "RATB1" => 0.0, "RATC1" => 0.0, "COD1" => 0,
    "CONT1" => 0, "RMA1" => 1.1, "RMI1" => 0.9,
    "VMA1" => 1.1, "VMI1" => 0.9, "NTP1" => 33,
    "TAB1" => 0, "CR1" => 0.0, "CX1" => 0.0, "CNXA1" => 0.0,
    "WINDV2" => nothing,
    "NOMV2" => 0.0, "ANG2" => 0.0, "RATA2" => 0.0,
    "RATB2" => 0.0, "RATC2" => 0.0, "COD2" => 0,
    "CONT2" => 0, "RMA2" => 1.1, "RMI2" => 0.9,
    "VMA2" => 1.1, "VMI2" => 0.9, "NTP2" => 33,
    "TAB2" => 0, "CR2" => 0.0, "CX2" => 0.0,
    "CNXA2" => 0.0,
    "WINDV3" => nothing,
    "NOMV3" => 0.0, "ANG3" => 0.0, "RATA3" => 0.0,
    "RATB3" => 0.0, "RATC3" => 0.0, "COD3" => 0,
    "CONT3" => 0, "RMA3" => 1.1, "RMI3" => 0.9,
    "VMA3" => 1.1, "VMI3" => 0.9, "NTP3" => 33,
    "TAB3" => 0, "CR3" => 0.0, "CX3" => 0.0,
    "CNXA3" => 0.0)

const _default_area_interchange = Dict("ISW" => 0, "PDES" => 0.0,
    "PTOL" => 10.0, "ARNAME" => "            ")

const _default_two_terminal_dc = Dict("MDC" => 0, "VCMOD" => 0.0,
    "RCOMP" => 0.0, "DELTI" => 0.0, "METER" => "I", "DCVMIN" => 0.0,
    "CCCITMX" => 20, "CCCACC" => 1.0, "TRR" => 1.0, "TAPR" => 1.0,
    "TMXR" => 1.5, "TMNR" => 0.51, "STPR" => 0.00625, "ICR" => 0, "IFR" => 0,
    "ITR" => 0, "IDR" => "1", "XCAPR" => 0.0, "TRI" => 1.0, "TAPI" => 1.0,
    "TMXI" => 1.5, "TMNI" => 0.51, "STPI" => 0.00625, "ICI" => 0, "IFI" => 0,
    "ITI" => 0, "IDI" => "1", "XCAPI" => 0.0)

const _default_vsc_dc = Dict("MDC" => 1, "O1" => nothing, "O2" => 0, "O3" => 0,
    "O4" => 0, "F1" => 1.0, "F2" => 1.0, "F3" => 1.0, "F4" => 1.0,
    "CONVERTER BUSES" => Dict("MODE" => 1, "ACSET" => 1.0, "ALOSS" => 1.0,
        "BLOSS" => 0.0, "MINLOSS" => 0.0, "SMAX" => 0.0, "IMAX" => 0.0,
        "PWF" => 1.0, "MAXQ" => 9999.0, "MINQ" => -9999.0, "REMOT" => 0,
        "RMPCT" => 100.0)
    )

const _default_impedance_correction = Dict("T1" => 0.0, "T2" => 0.0, "T3" => 0.0,
    "T4" => 0.0, "T5" => 0.0, "T6" => 0.0, "T7" => 0.0, "T8" => 0.0, "T9" => 0.0,
    "T10" => 0.0, "T11" => 0.0, "F1" => 0.0, "F2" => 0.0, "F3" => 0.0,
    "F4" => 0.0, "F5" => 0.0, "F6" => 0.0, "F7" => 0.0, "F8" => 0.0,
    "F9" => 0.0, "F10" => 0.0, "F11" => 0.0)

const _default_multi_term_dc = Dict("MDC" => 0, "VCMOD" => 0.0, "VCONVN" => 0,
    "CONV" => Dict("TR" => 1.0, "TAP" => 1.0, "TPMX" => 1.5, "TPMN" => 0.51,
        "TSTP" => 0.00625,"DCPF" => 1, "MARG" => 0.0, "CNVCOD" => 1),
    "DCBS" => Dict("IB" => 0.0, "AREA" => 1, "ZONE" => 1,
        "DCNAME" => "            ", "IDC2" => 0, "RGRND" => 0.0, "OWNER" => 1),
    "DCLN" => Dict("DCCKT" => 1, "MET" => 1, "LDC" => 0.0)
    )

const _default_multi_section = Dict("ID" => "&1", "MET" => 1)

const _default_zone = Dict("ZONAME" => "            ")

const _default_interarea = Dict("TRID" => 1, "PTRAN" => 0.0)

const _default_owner = Dict("OWNAME" => "            ")

const _default_facts = Dict("J" => 0, "MODE" => 1, "PDES" => 0.0, "QDES" => 0.0,
    "VSET" => 1.0, "SHMX" => 9999.0, "TRMX" => 9999.0, "VTMN" => 0.9,
    "VTMX" => 1.1, "VSMX" => 1.0, "IMX" => 0.0, "LINX" => 0.05,
    "RMPCT" => 100.0, "OWNER" => 1, "SET1" => 0.0, "SET2" => 0.0, "VSREF" => 0,
    "REMOT" => 0, "MNAME" => "")

const _default_switched_shunt = Dict("MODSW" => 1, "ADJM" => 0, "STAT" => 1,
    "VSWHI" => 1.0, "VSWLO" => 1.0, "SWREM" => 0, "RMPCT" => 100.0,
    "RMIDNT" => "", "BINIT" => 0.0,
    "N1" => 0.0, "N2" => 0.0, "N3" => 0.0, "N4" => 0.0, "N5" => 0.0,
    "N6" => 0.0, "N7" => 0.0, "N8" => 0.0,
    "B1" => 0.0, "B2" => 0.0, "B3" => 0.0, "B4" => 0.0, "B5" => 0.0,
    "B6" => 0.0, "B7" => 0.0, "B8" => 0.0)

const _default_gne_device = Dict("NTERM" => 1, "NREAL" => 0, "NINTG" => 0,
    "NCHAR" => 0, "STATUS" => 1, "OWNER" => nothing, "NMETR" => nothing,
    "REAL" => 0, "INTG" => nothing,
    "CHAR" => "1")

const _default_induction_machine = Dict("ID" => 1, "STAT" => 1, "SCODE" => 1,
    "DCODE" => 2, "AREA" => nothing, "ZONE" => nothing, "OWNER" => nothing,
    "TCODE" => 1, "BCODE" => 1, "MBASE" => nothing, "RATEKV" => 0.0,
    "PCODE" => 1, "H" => 1.0, "A" => 1.0, "B" => 1.0, "D" => 1.0, "E" => 1.0,
    "RA" => 0.0, "XA" => 0.0, "XM" => 2.5, "R1" => 999.0, "X1" => 999.0,
    "R2" => 999.0, "X2" => 999.0, "X3" => 0.0, "E1" => 1.0, "SE1" => 0.0,
    "E2" => 1.2, "SE2" => 0.0, "IA1" => 0.0, "IA2" => 0.0, "XAMULT" => 1)

const _pti_defaults = Dict("BUS" => _default_bus,
    "LOAD" => _default_load,
    "FIXED SHUNT" => _default_fixed_shunt,
    "GENERATOR" => _default_generator,
    "BRANCH" => _default_branch,
    "TRANSFORMER" => _default_transformer,
    "AREA INTERCHANGE" => _default_area_interchange,
    "TWO-TERMINAL DC" => _default_two_terminal_dc,
    "VOLTAGE SOURCE CONVERTER" => _default_vsc_dc,
    "IMPEDANCE CORRECTION" => _default_impedance_correction,
    "MULTI-TERMINAL DC" => _default_multi_term_dc,
    "MULTI-SECTION LINE" => _default_multi_section,
    "ZONE" => _default_zone,
    "INTER-AREA TRANSFER" => _default_interarea,
    "OWNER" => _default_owner,
    "FACTS CONTROL DEVICE" => _default_facts,
    "SWITCHED SHUNT" => _default_switched_shunt,
    "CASE IDENTIFICATION" => _default_case_identification,
    "GNE DEVICE" => _default_gne_device,
    "INDUCTION MACHINE" => _default_induction_machine
    )



function _correct_nothing_values!(data::Dict)
    if !haskey(data, "BUS")
        return
    end

    sbase = data["CASE IDENTIFICATION"][1]["SBASE"]
    bus_lookup = Dict(bus["I"] => bus for bus in data["BUS"])

    if haskey(data, "LOAD")
        for load in data["LOAD"]
            load_bus = bus_lookup[load["I"]]
            if load["AREA"] == nothing
                load["AREA"] = load_bus["AREA"]
            end
            if load["ZONE"] == nothing
                load["ZONE"] = load_bus["ZONE"]
            end
            if load["OWNER"] == nothing
                load["OWNER"] = load_bus["OWNER"]
            end
        end
    end

    if haskey(data, "GENERATOR")
        for gen in data["GENERATOR"]
            gen_bus = bus_lookup[gen["I"]]
            if haskey(gen, "OWNER") && gen["OWNER"] == nothing
                gen["OWNER"] = gen_bus["OWNER"]
            end
            if gen["MBASE"] == nothing
                gen["MBASE"] = sbase
            end
        end
    end

    if haskey(data, "BRANCH")
        for branch in data["BRANCH"]
            branch_bus = bus_lookup[branch["I"]]
            if haskey(branch, "OWNER") && branch["OWNER"] == nothing
                branch["OWNER"] = branch_bus["OWNER"]
            end
        end
    end

    if haskey(data, "TRANSFORMER")
        for transformer in data["TRANSFORMER"]
            transformer_bus = bus_lookup[transformer["I"]]
            for base_id in ["SBASE1-2", "SBASE2-3", "SBASE3-1"]
                if haskey(transformer, base_id) && transformer[base_id] == nothing
                    transformer[base_id] = sbase
                end
            end
            for winding_id in ["WINDV1", "WINDV2", "WINDV3"]
                if haskey(transformer, winding_id) &&  transformer[winding_id] == nothing
                    if transformer["CW"] == 2
                        transformer[winding_id] = transformer_bus["BASKV"]
                    else
                        transformer[winding_id] = 1.0
                    end
                end
            end
        end
    end

    #=
    if haskey(data, "VOLTAGE SOURCE CONVERTER")
        for mdc in data["VOLTAGE SOURCE CONVERTER"]
            # TODO 
            # "O1" => Expr(:call, :_get_component_property, data["BUS"], "OWNER", "I", get(get(component, "CONVERTER BUSES", [Dict()])[1], "IBUS", 0)),
        end
    end
    =#

    if haskey(data, "GNE DEVICE")
        for gne in data["GNE DEVICE"]
            gne_bus = bus_lookup[gne["I"]]
            if haskey(gne, "OWNER") && gne["OWNER"] == nothing
                gne["OWNER"] = gne_bus["OWNER"]
            end
            if haskey(gne, "NMETR") && gne["NMETR"] == nothing
                gne["NMETR"] = gne_bus["NTERM"]
            end
        end
    end

    if haskey(data, "INDUCTION MACHINE")
        for indm in data["INDUCTION MACHINE"]
            indm_bus = bus_lookup[indm["I"]]
            if indm["AREA"] == nothing
                indm["AREA"] = indm_bus["AREA"]
            end
            if indm["ZONE"] == nothing
                indm["ZONE"] = indm_bus["ZONE"]
            end
            if indm["OWNER"] == nothing
                indm["OWNER"] = indm_bus["OWNER"]
            end
            if indm["MBASE"] == nothing
                indm["MBASE"] = sbase
            end
        end
    end
end



"experimental functional method for parsing elements and setting defaults"
function _parse_elements(elements::Array, dtypes::Array, defaults::Dict, section::AbstractString)
    data = Dict{String,Any}()

    if length(elements) > length(dtypes)
        Memento.warn(_LOGGER, "ignoring $(length(elements) - length(dtypes)) extra values in section $section, only $(length(dtypes)) items are defined")
        elements = elements[1:length(dtypes)]
    end

    for (i,element) in enumerate(elements)
        field, dtype = dtypes[i]

        element = strip(element)

        if dtype == String
            if startswith(element, "'") && endswith(element, "'")
                data[field] = element[2:end-1]
            else
                data[field] = element
            end
        else
            if length(element) <= 0
                # this will be set to a default in the cleanup phase
                data[field] = nothing
            else
                try
                    data[field] = parse(dtype, element)
                catch message
                    if isa(message, Meta.ParseError)
                        data[field] = element
                    else
                        Memento.error(_LOGGER, "value '$element' for $field in section $section is not of type $dtype.")
                    end
                end
            end
        end
    end

    if length(elements) < length(dtypes)
        for (field, dtype) in dtypes[length(elements):end]
            data[field] = defaults[field]
            #=
            if length(missing_fields) > 0
                for field in missing_fields
                    data[field] = ""
                end
                missing_str = join(missing_fields, ", ")
                if !(section == "SWITCHED SHUNT" && startswith(missing_str, "N")) &&
                    !(section == "MULTI-SECTION LINE" && startswith(missing_str, "DUM")) &&
                    !(section == "IMPEDANCE CORRECTION" && startswith(missing_str, "T"))
                    Memento.warn(_LOGGER, "The following fields in $section are missing: $missing_str")
                end
            end
            =#
        end
    end

    return data
end


"""
    _parse_line_element!(data, elements, section)

Internal function. Parses a single "line" of data elements from a PTI file, as
given by `elements` which is an array of the line, typically split at `,`.
Elements are parsed into data types given by `section` and saved into `data::Dict`.
"""
function _parse_line_element!(data::Dict, elements::Array, section::AbstractString)
    missing_fields = []
    for (i, (field, dtype)) in enumerate(_pti_dtypes[section])
        if i > length(elements)
            Memento.debug(_LOGGER, "Have run out of elements in $section at $field")
            push!(missing_fields, field)
            continue
        else
            element = strip(elements[i])
        end

        try
            if dtype != String && element != ""
                data[field] = parse(dtype, element)
            else
                if dtype == String && startswith(element, "'") && endswith(element, "'")
                    data[field] = chop(element[nextind(element,1):end])
                else
                    data[field] = element
                end
            end
        catch message
            if isa(message, Meta.ParseError)
                data[field] = element
            else
                Memento.error(_LOGGER, "value '$element' for $field in section $section is not of type $dtype.")
            end
        end
    end

    if length(missing_fields) > 0
        for field in missing_fields
            data[field] = ""
        end
        missing_str = join(missing_fields, ", ")
        if !(section == "SWITCHED SHUNT" && startswith(missing_str, "N")) &&
            !(section == "MULTI-SECTION LINE" && startswith(missing_str, "DUM")) &&
            !(section == "IMPEDANCE CORRECTION" && startswith(missing_str, "T"))
            Memento.warn(_LOGGER, "The following fields in $section are missing: $missing_str")
        end
    end
end



const _comment_split = r"(?!\B[\'][^\']*)[\/](?![^\']*[\']\B)"
const _split_string = r",(?=(?:[^']*'[^']*')*[^']*$)"

"""
    _get_line_elements(line)

Internal function. Uses regular expressions to extract all separate data
elements from a line of a PTI file and populate them into an `Array{String}`.
Comments, typically indicated at the end of a line with a `'/'` character,
are also extracted separately, and `Array{Array{String}, String}` is returned.
"""
function _get_line_elements(line::AbstractString)
    if count(i->(i=="'"), line) % 2 == 1
        throw(Memento.error(_LOGGER, "There are an uneven number of single-quotes in \"{line}\", the line cannot be parsed."))
    end

    line_comment = split(line, _comment_split, limit=2)
    line = strip(line_comment[1])
    comment = length(line_comment) > 1 ? strip(line_comment[2]) : ""

    #elements = [strip(e) for e in split(line, _split_string)]
    elements = split(line, _split_string)

    # results in a 20% memory overhead
    #Memento.debug(_LOGGER, "$line")
    #Memento.debug(_LOGGER, "$comment")
    #Memento.debug(_LOGGER, "$elements")

    return (elements, comment)
end


"""
    _parse_pti_data(data_string, sections)

Internal function. Parse a PTI raw file into a `Dict`, given the
`data_string` of the file and a list of the `sections` in the PTI
file (typically given by default by `get_pti_sections()`.
"""
function _parse_pti_data(data_io::IO)
    sections = deepcopy(_pti_sections)
    data_lines = readlines(data_io)
    skip_lines = 0
    skip_sublines = 0
    subsection = ""

    pti_data = Dict{String,Array{Dict}}()

    section = popfirst!(sections)
    section_data = Dict{String,Any}()

    for (line_number, line) in enumerate(data_lines)
        #Memento.debug(_LOGGER, "$line_number: $line")

        (elements, comment) = _get_line_elements(line)

        first_element = strip(elements[1])
        if line_number > 3 && length(elements) != 0 && first_element == "Q"
            break
        elseif line_number > 3 && length(elements) != 0 && first_element == "0"
            if line_number == 4
                section = popfirst!(sections)
            end

            if length(elements) > 1
                Memento.warn(_LOGGER, "At line $line_number, new section started with '0', but additional non-comment data is present. Pattern '^\\s*0\\s*[/]*.*' is reserved for section start/end.")
            elseif length(comment) > 0
                Memento.debug(_LOGGER, "At line $line_number, switched to $section")
            end

            if !isempty(sections)
                section = popfirst!(sections)
            end

            continue
        else
            if line_number == 4
                section = popfirst!(sections)
                section_data = Dict{String,Any}()
            end

            if skip_lines > 0
                skip_lines -= 1
                continue
            end

            Memento.debug(_LOGGER, join(["Section:", section], " "))
            if !(section in ["CASE IDENTIFICATION","TRANSFORMER","VOLTAGE SOURCE CONVERTER","MULTI-TERMINAL DC","TWO-TERMINAL DC","GNE DEVICE"])
                section_data = Dict{String,Any}()

                try
                    #dtypes = _pti_dtypes[section]
                    #defaults = _pti_defaults[section]
                    #section_data_tmp = _parse_elements(elements, dtypes, defaults, section)
                    _parse_line_element!(section_data, elements, section)
                catch message
                    throw(Memento.error(_LOGGER, "Parsing failed at line $line_number: $(sprint(showerror, message))"))
                end

            elseif section == "CASE IDENTIFICATION"
                if line_number == 1
                    try
                        _parse_line_element!(section_data, elements, section)
                    catch message
                        throw(Memento.error(_LOGGER, "Parsing failed at line $line_number: $(sprint(showerror, message))"))
                    end

                    if section_data["REV"] != "" && section_data["REV"] < 33
                        Memento.warn(_LOGGER, "Version $(section_data["REV"]) of PTI format is unsupported, parser may not function correctly.")
                    end
                else
                    section_data["Comment_Line_$(line_number - 1)"] = line
                end

                if line_number < 3
                    continue
                end

            elseif section == "TRANSFORMER"
                section_data = Dict{String,Any}()
                if parse(Int64, _get_line_elements(line)[1][3]) == 0 # two winding transformer
                    winding = "TWO-WINDING"
                    skip_lines = 3
                elseif parse(Int64, _get_line_elements(line)[1][3]) != 0 # three winding transformer
                    winding = "THREE-WINDING"
                    skip_lines = 4
                else
                    Memento.error(_LOGGER, "Cannot detect type of Transformer")
                end

                try
                    for transformer_line in 0:4
                        if transformer_line == 0
                            temp_section = section
                        else
                            temp_section = join([section, winding, "LINE", transformer_line], " ")
                        end

                        if winding == "TWO-WINDING" && transformer_line == 4
                            break
                        else
                            elements = _get_line_elements(data_lines[line_number + transformer_line])[1]
                            _parse_line_element!(section_data, elements, temp_section)
                        end
                    end
                catch message
                    throw(Memento.error(_LOGGER, "Parsing failed at line $line_number: $(sprint(showerror, message))"))
                end

            elseif section == "VOLTAGE SOURCE CONVERTER"
                if length(_get_line_elements(line)[1]) == 11
                    section_data = Dict{String,Any}()
                    try
                        _parse_line_element!(section_data, elements, section)
                    catch message
                        throw(Memento.error(_LOGGER, "Parsing failed at line $line_number: $(sprint(showerror, message))"))
                    end
                    skip_sublines = 2
                    continue

                elseif skip_sublines > 0
                    skip_sublines -= 1
                    subsection_data = Dict{String,Any}()

                    for (field, dtype) in _pti_dtypes["$section SUBLINES"]
                        element = popfirst!(elements)
                        if element != ""
                            subsection_data[field] = parse(dtype, element)
                        else
                            subsection_data[field] = ""
                        end
                    end

                    if haskey(section_data, "CONVERTER BUSES")
                        push!(section_data["CONVERTER BUSES"], subsection_data)
                    else
                        section_data["CONVERTER BUSES"] = [subsection_data]
                        continue
                    end
                end

            elseif section == "TWO-TERMINAL DC"
                section_data = Dict{String,Any}()
                if length(_get_line_elements(line)[1]) == 12
                    (elements, comment) = _get_line_elements(join(data_lines[line_number:line_number + 2], ','))
                    skip_lines = 2
                end

                try
                    _parse_line_element!(section_data, elements, section)
                catch message
                    throw(Memento.error(_LOGGER, "Parsing failed at line $line_number: $(sprint(showerror, message))"))
                end

            elseif section == "MULTI-TERMINAL DC"
                if skip_sublines == 0
                    section_data = Dict{String,Any}()
                    try
                        _parse_line_element!(section_data, elements, section)
                    catch message
                        throw(Memento.error(_LOGGER, "Parsing failed at line $line_number: $(sprint(showerror, message))"))
                    end

                    if section_data["NCONV"] > 0
                        skip_sublines = section_data["NCONV"]
                        subsection = "NCONV"
                        continue
                    elseif section_data["NDCBS"] > 0
                        skip_sublines = section_data["NDCBS"]
                        subsection = "NDCBS"
                        continue
                    elseif section_data["NDCLN"] > 0
                        skip_sublines = section_data["NDCLN"]
                        subsection = "NDCLN"
                        continue
                    end
                end

                if skip_sublines > 0
                    skip_sublines -= 1

                    subsection_data = Dict{String,Any}()
                    try
                        _parse_line_element!(subsection_data, elements, "$section $subsection")
                    catch message
                        throw(Memento.error(_LOGGER, "Parsing failed at line $line_number: $(sprint(showerror, message))"))
                    end

                    if haskey(section_data, "$(subsection[2:end])")
                        section_data["$(subsection[2:end])"] = push!(section_data["$(subsection[2:end])"], subsection_data)
                        if skip_sublines > 0 && subsection != "NDCLN"
                            continue
                        end
                    else
                        section_data["$(subsection[2:end])"] = [subsection_data]
                        if skip_sublines > 0 && subsection != "NDCLN"
                            continue
                        end
                    end

                    if skip_sublines == 0 && subsection != "NDCLN"
                        if subsection == "NDCBS"
                            skip_sublines = section_data["NDCLN"]
                            subsection = "NDCLN"
                            continue
                        elseif subsection == "NCONV"
                            skip_sublines = section_data["NDCBS"]
                            subsection = "NDCBS"
                            continue
                        end
                    elseif skip_sublines == 0 && subsection == "NDCLN"
                        subsection = ""
                    else
                        continue
                    end
                end

            elseif section == "GNE DEVICE"
                # TODO: handle multiple lines of GNE Device
                Memento.warn(_LOGGER, "GNE DEVICE parsing is not supported.")
            end
        end
        if subsection != ""
            Memento.debug(_LOGGER, "appending data")
        end

        if haskey(pti_data, section)
            push!(pti_data[section], section_data)
        else
            pti_data[section] = [section_data]
        end
    end

    _populate_defaults!(pti_data)
    _correct_nothing_values!(pti_data)

    return pti_data
end


"""
    parse_pti(filename::String)

Open PTI raw file given by `filename`, returning a `Dict` of the data parsed
into the proper types.
"""
function parse_pti(filename::String)::Dict
    pti_data = open(filename) do f
        parse_pti(f)
    end

    return pti_data
end


"""
    parse_pti(io::IO)

Reads PTI data in `io::IO`, returning a `Dict` of the data parsed into the
proper types.
"""
function parse_pti(io::IO)::Dict
    pti_data = _parse_pti_data(io)
    try
        pti_data["CASE IDENTIFICATION"][1]["NAME"] = match(r"^\<file\s[\/\\]*(?:.*[\/\\])*(.*)\.raw\>$", lowercase(io.name)).captures[1]
    catch
        throw(Memento.error(_LOGGER, "This file is unrecognized and cannot be parsed"))
    end

    return pti_data
end



"""
    _populate_defaults!(pti_data)

Internal function. Populates empty fields with PSS(R)E PTI v33 default values
"""
function _populate_defaults!(data::Dict)
    for section in _pti_sections
        if haskey(data, section)
            component_defaults = _pti_defaults[section]
            for component in data[section]
                for (field, field_value) in component
                    if isa(field_value, Array)
                        sub_component_defaults = component_defaults[field]
                        for sub_component in field_value
                            for (sub_field, sub_field_value) in sub_component
                                if sub_field_value == ""
                                    try
                                        sub_component[sub_field] = sub_component_defaults[sub_field]
                                    catch msg
                                        if isa(msg, KeyError)
                                            Memento.warn(_LOGGER, "'$sub_field' in '$field' in '$section' has no default value")
                                        else
                                            rethrow(msg)
                                        end
                                    end
                                end
                            end
                        end
                    elseif field_value == "" && !(field in ["Comment_Line_1", "Comment_Line_2"]) && !startswith(field, "DUM")
                        try
                            component[field] = component_defaults[field]
                        catch msg
                            if isa(msg, KeyError)
                                Memento.warn(_LOGGER, "'$field' in '$section' has no default value")
                            else
                                rethrow(msg)
                            end
                        end
                    end
                end
            end
        end
    end
end


