##### GOC Data Tools #####

# shared across C1 and C2


##### Contingency Description Data File (.con) #####

# OPEN BRANCH FROM BUS *I TO BUS *J CIRCUIT *1CKT
const _branch_contigency_structure = [
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
const _generator_contigency_structure = [
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
        #line_tokens = split(strip(_remove_psse_comment(line)))
        line_tokens = split(strip(line))
        #println(line_tokens)
        append!(tokens, line_tokens)
    end

    #println(tokens)

    token_idx = 1
    while token_idx <= length(tokens)
        token = tokens[token_idx]
        if token == "END"
            debug(_LOGGER, "end of contingency file found")
            break
        elseif token == "CONTINGENCY"
            # start reading contingencies

            contingency_name = tokens[token_idx+1]
            debug(_LOGGER, "reading contingency $(contingency_name)")

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
                    error(_LOGGER, "incorrect branch contingency structure: $(branch_tokens)")
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
                    error(_LOGGER, "incorrect generator contingency structure: $(generator_tokens)")
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
                warn(_LOGGER, "no action provided for contingency $(contingency_name)")
                token_idx -= 1
            else
                warn(_LOGGER, "unrecognized token $(token)")
            end

            token_idx += 1
            token = tokens[token_idx]
            if token != "END"
                error(_LOGGER, "expected END token at end of CONTINGENCY, got $(token)")
            end
        else
            warn(_LOGGER, "unrecognized token $(token)")
        end
        token_idx += 1
    end

    return con_lists
end
