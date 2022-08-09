const C1_PG_LOSS_TOL = 1e-6

function compute_susceptance_matrix_inv_size(network)
    b0 = Base.gc_bytes()
    b_inv = _PM.calc_susceptance_matrix_inv(network)
    return Base.gc_bytes() - b0
end


# a global network variable used for iterative computations
c1_network_global = Dict{String,Any}()

# a global contingency list used for iterative computations
c1_contingency_order_global = []

""
function load_c1_network_global(con_file, inl_file, raw_file, rop_file, scenario_id)
    info(_LOGGER, "skipping goc and power models data warnings")
    pm_logger_level = getlevel(getlogger(PowerModels))
    goc_logger_level = getlevel(_LOGGER)

    setlevel!(getlogger(PowerModels), "error")
    setlevel!(_LOGGER, "error")

    goc_data = parse_c1_files(con_file, inl_file, raw_file, rop_file, scenario_id=scenario_id)
    global c1_network_global = build_c1_pm_model(goc_data)
    global c1_contingency_order_global = contingency_order(c1_network_global)

    setlevel!(getlogger(PowerModels), pm_logger_level)
    setlevel!(_LOGGER, goc_logger_level)

    return 0
end


function write_c1_contingencies(network; output_dir="", file_name="contingencies.txt")
    if length(output_dir) > 0
        file_path = joinpath(output_dir, file_name)
    else
        file_path = file_name
    end

    open(file_path, "w") do cont_file
        for c in network["gen_contingencies_active"]
            write(cont_file, "$(c.label)\n")
        end
        for c in network["branch_contingencies_active"]
            write(cont_file, "$(c.label)\n")
        end
    end
end

function write_c1_active_flow_cuts(network; output_dir="", file_name="active_flow_cuts.txt")
    if length(output_dir) > 0
        file_path = joinpath(output_dir, file_name)
    else
        file_path = file_name
    end

    open(file_path, "w") do cont_file
        for cut in network["gen_flow_cuts"]
            write(cont_file, "$(cut.cont_label), gen, $(cut.branch_id), $(cut.rating_level)\n")
        end

        for cut in network["branch_flow_cuts"]
            write(cont_file, "$(cut.cont_label), branch, $(cut.branch_id), $(cut.rating_level)\n")
        end
    end
end

function read_c1_active_flow_cuts(;output_dir="", file_name="active_flow_cuts.txt")
    if length(output_dir) > 0
        file_path = joinpath(output_dir, file_name)
    else
        file_path = file_name
    end

    if isfile(file_path)
        info(_LOGGER, "loading flow cuts file: $(file_path)")
        return parse_c1_flow_cuts_file(file_path)
    else
        info(_LOGGER, "flow cuts file not found: $(file_path)")
        return []
    end
end

function parse_c1_flow_cuts_file(file::String)
    open(file) do io
        return parse_c1_flow_cuts_file(io)
    end
end

function parse_c1_flow_cuts_file(io::IO)
    cuts_list = []

    for line in readlines(io)
        if length(strip(line)) == 0
            warn(_LOGGER, "skipping blank line in cuts file")
            continue
        end
        line_parts = split(line, ",")
        if length(line_parts) != 4
            warn(_LOGGER, "skipping ill formated line\n   $(line)")
            continue
        end

        cont_label = strip(line_parts[1])
        cont_type = strip(line_parts[2])
        branch_id = parse(Int, line_parts[3])
        rating_level = parse(Float64, line_parts[4])

        push!(cuts_list, (cont_label=cont_label, cont_type=cont_type, branch_id=branch_id, rating_level=rating_level))
    end

    return cuts_list
end



"ranks generator contingencies and down selects based on evaluation limits"
function calc_c1_gen_contingency_subset(network::Dict{String,<:Any}; gen_eval_limit=length(network["gen_contingencies"]))
    gen_cap = Dict(gen["index"] => sqrt(max(abs(gen["pmin"]), abs(gen["pmax"]))^2 + max(abs(gen["qmin"]), abs(gen["qmax"]))^2) for (i,gen) in network["gen"])
    gen_contingencies = sort(network["gen_contingencies"], rev=true, by=x -> gen_cap[x.idx])

    gen_cont_limit = min(gen_eval_limit, length(network["gen_contingencies"]))
    gen_contingencies = gen_contingencies[1:gen_cont_limit]

    return gen_contingencies
end

"ranks branch contingencies and down selects based on evaluation limits"
function calc_c1_branch_contingency_subset(network::Dict{String,<:Any}; branch_eval_limit=length(network["branch_contingencies"]))
    line_imp_mag = Dict(branch["index"] => branch["rate_a"]*sqrt(branch["br_r"]^2 + branch["br_x"]^2) for (i,branch) in network["branch"])
    branch_contingencies = sort(network["branch_contingencies"], rev=true, by=x -> line_imp_mag[x.idx])

    branch_cont_limit = min(branch_eval_limit, length(network["branch_contingencies"]))
    branch_contingencies = branch_contingencies[1:branch_cont_limit]

    return branch_contingencies
end


function calc_c1_branch_ptdf_single(am::_PM.AdmittanceMatrix, ref_bus::Int, branch::Dict{String,<:Any})
    branch_ptdf = Dict{Int,Any}()
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]

    b = imag(inv(branch["br_r"] + im * branch["br_x"]))

    va_fr = _PM.injection_factors_va(am, ref_bus, f_bus)
    va_to = _PM.injection_factors_va(am, ref_bus, t_bus)

    # convert bus injection functions to PTDF style
    bus_injection = Dict(i => -b*(get(va_fr, i, 0.0) - get(va_to, i, 0.0)) for i in union(keys(va_fr), keys(va_to)))

    return bus_injection
end





""
function check_c1_contingency_violations_remote(cont_range, output_dir, contingency_limit=1, solution_file="solution1.txt")
    if length(c1_network_global) <= 0 || length(c1_contingency_order_global) <= 0
        error(_LOGGER, "check_contingencies_branch_flow_remote called before load_c1_network_global")
    end

    sol = read_c1_solution1(c1_network_global, output_dir=output_dir, state_file=solution_file)
    _PM.update_data!(c1_network_global, sol)

    network = copy(c1_network_global)
    contingencies = c1_contingency_order_global[cont_range]
    network["gen_contingencies"] = [c for c in contingencies if c.type == "gen"]
    network["branch_contingencies"] = [c for c in contingencies if c.type == "branch"]

    contingencies = check_c1_contingency_violations(network, contingency_limit=contingency_limit)

    return contingencies
end


"""
Checks a given operating point against the contingencies to look for branch
flow violations.  The DC Power Flow approximation is used for flow simulation.
Returns a list of contingencies where a violation is found.
"""
function check_c1_contingency_violations(network;
        gen_contingency_limit=10, branch_contingency_limit=10, contingency_limit=typemax(Int64),
        gen_eval_limit=typemax(Int64), branch_eval_limit=typemax(Int64), sm_threshold=0.01)

    if _IM.ismultinetwork(network)
        error(_LOGGER, "the branch flow cut generator can only be used on single networks")
    end
    time_contingencies_start = time()


    network_lal = deepcopy(network) #lal -> losses as loads

    ref_bus_id = _PM.reference_bus(network_lal)["index"]

    gen_pg_init = Dict(i => gen["pg"] for (i,gen) in network_lal["gen"])

    load_active = Dict(i => load for (i,load) in network_lal["load"] if load["status"] != 0)

    pd_total = sum(load["pd"] for (i,load) in load_active)
    p_losses = sum(gen["pg"] for (i,gen) in network_lal["gen"] if gen["gen_status"] != 0) - pd_total
    p_delta = 0.0

    if p_losses > C1_PG_LOSS_TOL
        load_count = length(load_active)
        p_delta = p_losses/load_count
        for (i,load) in load_active
            load["pd"] += p_delta
        end
        warn(_LOGGER, "active power losses found $(p_losses) increasing loads by $(p_delta)")
    end



    gen_contingencies = calc_c1_gen_contingency_subset(network_lal, gen_eval_limit=gen_eval_limit)
    branch_contingencies = calc_c1_branch_contingency_subset(network_lal, branch_eval_limit=branch_eval_limit)

    gen_cuts = []
    for (i,cont) in enumerate(gen_contingencies)
        if length(gen_cuts) >= gen_contingency_limit
            info(_LOGGER, "hit gen flow cut limit $(gen_contingency_limit)")
            break
        end
        if length(gen_cuts) >= contingency_limit
            info(_LOGGER, "hit total cut limit $(contingency_limit)")
            break
        end
        #info(_LOGGER, "working on ($(i)/$(gen_eval_limit)/$(gen_cont_total)): $(cont.label)")

        for (i,gen) in network_lal["gen"]
            gen["pg"] = gen_pg_init[i]
        end

        cont_gen = network_lal["gen"]["$(cont.idx)"]
        pg_lost = cont_gen["pg"]

        cont_gen["gen_status"] = 0
        cont_gen["pg"] = 0.0


        gen_bus = network_lal["bus"]["$(cont_gen["gen_bus"])"]
        gen_set = network_lal["area_gens"][gen_bus["area"]]

        gen_active = Dict(i => gen for (i,gen) in network_lal["gen"] if gen["index"] != cont.idx && gen["index"] in gen_set && gen["gen_status"] != 0)

        alpha_gens = [gen["alpha"] for (i,gen) in gen_active]
        if length(alpha_gens) == 0 || isapprox(sum(alpha_gens), 0.0)
            warn(_LOGGER, "no available active power response in cont $(cont.label), active gens $(length(alpha_gens))")
            continue
        end

        alpha_total = sum(alpha_gens)
        delta = pg_lost/alpha_total
        network_lal["delta"] = delta
        #info(_LOGGER, "$(pg_lost) - $(alpha_total) - $(delta)")

        for (i,gen) in gen_active
            gen["pg"] += gen["alpha"]*delta
        end

        try
            solution = _PM.compute_dc_pf(network_lal)["solution"]
            _PM.update_data!(network_lal, solution)
        catch exception
            warn(_LOGGER, "linear solve failed on $(cont.label)")
            continue
        end

        flow = _PM.calc_branch_flow_dc(network_lal)
        _PM.update_data!(network_lal, flow)


        vio = calc_c1_violations(network_lal, network_lal)

        #info(_LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            info(_LOGGER, "adding contingency $(cont.label) due to constraint flow violations $(vio.sm)")
            push!(gen_cuts, cont)
        end

        cont_gen["gen_status"] = 1
        cont_gen["pg"] = pg_lost
        network_lal["delta"] = 0.0
    end


    branch_cuts = []
    for (i,cont) in enumerate(branch_contingencies)
        if length(branch_cuts) >= branch_contingency_limit
            info(_LOGGER, "hit branch flow cut limit $(branch_contingency_limit)")
            break
        end
        if length(gen_cuts) + length(branch_cuts) >= contingency_limit
            info(_LOGGER, "hit total cut limit $(contingency_limit)")
            break
        end

        #info(_LOGGER, "working on ($(i)/$(branch_eval_limit)/$(branch_cont_total)): $(cont.label)")

        cont_branch = network_lal["branch"]["$(cont.idx)"]
        cont_branch["br_status"] = 0

        try
            solution = _PM.compute_dc_pf(network_lal)["solution"]
            _PM.update_data!(network_lal, solution)
        catch exception
            warn(_LOGGER, "linear solve failed on $(cont.label)")
            cont_branch["br_status"] = 1
            continue
        end

        flow = _PM.calc_branch_flow_dc(network_lal)
        _PM.update_data!(network_lal, flow)

        vio = calc_c1_violations(network_lal, network_lal)

        #info(_LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            info(_LOGGER, "adding contingency $(cont.label) due to constraint flow violations $(vio.sm)")
            push!(branch_cuts, cont)
        end

        cont_branch["br_status"] = 1
    end


    if p_delta != 0.0
        warn(_LOGGER, "re-adjusting loads by $(-p_delta)")
        for (i,load) in load_active
            load["pd"] -= p_delta
        end
    end

    time_contingencies = time() - time_contingencies_start
    info(_LOGGER, "contingency eval time: $(time_contingencies)")

    return (gen_contingencies=gen_cuts, branch_contingencies=branch_cuts)
end





""
function check_c1_contingencies_branch_power_remote(cont_range, output_dir, cut_limit=1, solution_file="solution1.txt")
    if length(c1_network_global) <= 0 || length(c1_contingency_order_global) <= 0
        error(_LOGGER, "check_contingencies_branch_flow_remote called before load_c1_network_global")
    end

    sol = read_c1_solution1(c1_network_global, output_dir=output_dir, state_file=solution_file)
    _PM.update_data!(c1_network_global, sol)

    active_cuts = read_c1_active_flow_cuts(output_dir=output_dir)
    gen_flow_cuts = []
    branch_flow_cuts = []
    for cut in active_cuts
        if cut.cont_type == "gen"
            push!(gen_flow_cuts, cut)
        elseif cut.cont_type == "branch"
            push!(branch_flow_cuts, cut)
        else
            warn(_LOGGER, "unknown contingency type in cut $(cut)")
        end
    end

    network = copy(c1_network_global)
    contingencies = c1_contingency_order_global[cont_range]
    network["gen_contingencies"] = [c for c in contingencies if c.type == "gen"]
    network["branch_contingencies"] = [c for c in contingencies if c.type == "branch"]

    cuts = check_c1_contingencies_branch_power(network, total_cut_limit=cut_limit, gen_flow_cuts=gen_flow_cuts, branch_flow_cuts=branch_flow_cuts)

    return cuts
end


"""
Checks a given operating point against the contingencies to look for branch
flow violations.  The DC Power Flow approximation is used for flow simulation.
If a violation is found, computes a PTDF cut based on bus injections.  Uses the
participation factor based generator response model from the ARPA-e GOC
Challenge 1 specification.
"""
function check_c1_contingencies_branch_power(network;
        gen_flow_cut_limit=10, branch_flow_cut_limit=10, total_cut_limit=typemax(Int64),
        gen_eval_limit=typemax(Int64), branch_eval_limit=typemax(Int64), sm_threshold=0.01,
        gen_flow_cuts=[], branch_flow_cuts=[]
        )

    if _IM.ismultinetwork(network)
        error(_LOGGER, "the branch flow cut generator can only be used on single networks")
    end
    time_contingencies_start = time()

    gen_cuts_active = Dict()
    for gen_cut in gen_flow_cuts
        if !haskey(gen_cuts_active, gen_cut.cont_label)
            gen_cuts_active[gen_cut.cont_label] = Set{Int}()
        end
        push!(gen_cuts_active[gen_cut.cont_label], gen_cut.branch_id)
    end

    branch_cuts_active = Dict()
    for branch_cut in branch_flow_cuts
        if !haskey(branch_cuts_active, branch_cut.cont_label)
            branch_cuts_active[branch_cut.cont_label] = Set{Int}()
        end
        push!(branch_cuts_active[branch_cut.cont_label], branch_cut.branch_id)
    end


    network_lal = deepcopy(network) #lal -> losses as loads

    ref_bus_id = _PM.reference_bus(network_lal)["index"]

    gen_pg_init = Dict(i => gen["pg"] for (i,gen) in network_lal["gen"])

    load_active = Dict(i => load for (i,load) in network_lal["load"] if load["status"] != 0)

    pd_total = sum(load["pd"] for (i,load) in load_active)
    p_losses = sum(gen["pg"] for (i,gen) in network_lal["gen"] if gen["gen_status"] != 0) - pd_total
    p_delta = 0.0

    if p_losses > C1_PG_LOSS_TOL
        load_count = length(load_active)
        p_delta = p_losses/load_count
        for (i,load) in load_active
            load["pd"] += p_delta
        end
        warn(_LOGGER, "active power losses found $(p_losses) increasing loads by $(p_delta)")
    end


    gen_contingencies = calc_c1_gen_contingency_subset(network_lal, gen_eval_limit=gen_eval_limit)
    branch_contingencies = calc_c1_branch_contingency_subset(network_lal, branch_eval_limit=branch_eval_limit)

    gen_cuts = []
    for (i,cont) in enumerate(gen_contingencies)
        if length(gen_cuts) >= gen_flow_cut_limit
            info(_LOGGER, "hit gen flow cut limit $(gen_flow_cut_limit)")
            break
        end
        if length(gen_cuts) >= total_cut_limit
            info(_LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end
        #info(_LOGGER, "working on ($(i)/$(gen_eval_limit)/$(gen_cont_total)): $(cont.label)")

        for (i,gen) in network_lal["gen"]
            gen["pg"] = gen_pg_init[i]
        end

        cont_gen = network_lal["gen"]["$(cont.idx)"]
        pg_lost = cont_gen["pg"]

        cont_gen["gen_status"] = 0
        cont_gen["pg"] = 0.0


        gen_bus = network_lal["bus"]["$(cont_gen["gen_bus"])"]
        gen_set = network_lal["area_gens"][gen_bus["area"]]

        gen_active = Dict(i => gen for (i,gen) in network_lal["gen"] if gen["index"] != cont.idx && gen["index"] in gen_set && gen["gen_status"] != 0)

        alpha_gens = [gen["alpha"] for (i,gen) in gen_active]
        if length(alpha_gens) == 0 || isapprox(sum(alpha_gens), 0.0)
            warn(_LOGGER, "no available active power response in cont $(cont.label), active gens $(length(alpha_gens))")
            continue
        end

        alpha_total = sum(alpha_gens)
        delta = pg_lost/alpha_total
        network_lal["delta"] = delta
        #info(_LOGGER, "$(pg_lost) - $(alpha_total) - $(delta)")

        for (i,gen) in gen_active
            gen["pg"] += gen["alpha"]*delta
        end

        try
            solution = _PM.compute_dc_pf(network_lal)["solution"]
            _PM.update_data!(network_lal, solution)
        catch exception
            warn(_LOGGER, "linear solve failed on $(cont.label)")
            continue
        end

        flow = _PM.calc_branch_flow_dc(network_lal)
        _PM.update_data!(network_lal, flow)


        vio = calc_c1_violations(network_lal, network_lal)

        #info(_LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vios = branch_c1_violations_sorted(network_lal, network_lal)
            branch_vio = branch_vios[1]

            if !haskey(gen_cuts_active, cont.label) || !(branch_vio.branch_id in gen_cuts_active[cont.label])
                info(_LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                am = _PM.calc_susceptance_matrix(network_lal)
                branch = network_lal["branch"]["$(branch_vio.branch_id)"]

                bus_injection = calc_c1_branch_ptdf_single(am, ref_bus_id, branch)
                cut = (gen_id=cont.idx, cont_label=cont.label, branch_id=branch_vio.branch_id, rating_level=1.0, bus_injection=bus_injection)
                push!(gen_cuts, cut)
            else
                warn(_LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end

        end

        cont_gen["gen_status"] = 1
        cont_gen["pg"] = pg_lost
        network_lal["delta"] = 0.0
    end


    branch_cuts = []
    for (i,cont) in enumerate(branch_contingencies)
        if length(branch_cuts) >= branch_flow_cut_limit
            info(_LOGGER, "hit branch flow cut limit $(branch_flow_cut_limit)")
            break
        end
        if length(gen_cuts) + length(branch_cuts) >= total_cut_limit
            info(_LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end

        #info(_LOGGER, "working on ($(i)/$(branch_eval_limit)/$(branch_cont_total)): $(cont.label)")

        cont_branch = network_lal["branch"]["$(cont.idx)"]
        cont_branch["br_status"] = 0

        try
            solution = _PM.compute_dc_pf(network_lal)["solution"]
            _PM.update_data!(network_lal, solution)
        catch exception
            warn(_LOGGER, "linear solve failed on $(cont.label)")
            cont_branch["br_status"] = 1
            continue
        end

        flow = _PM.calc_branch_flow_dc(network_lal)
        _PM.update_data!(network_lal, flow)

        vio = calc_c1_violations(network_lal, network_lal)

        #info(_LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vio = branch_c1_violations_sorted(network_lal, network_lal)[1]
            if !haskey(branch_cuts_active, cont.label) || !(branch_vio.branch_id in branch_cuts_active[cont.label])
                info(_LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                am = _PM.calc_susceptance_matrix(network_lal)
                branch = network_lal["branch"]["$(branch_vio.branch_id)"]

                bus_injection = calc_c1_branch_ptdf_single(am, ref_bus_id, branch)
                cut = (cont_label=cont.label, branch_id=branch_vio.branch_id, rating_level=1.0, bus_injection=bus_injection)
                push!(branch_cuts, cut)
            else
                warn(_LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end
        end

        cont_branch["br_status"] = 1
    end


    if p_delta != 0.0
        warn(_LOGGER, "re-adjusting loads by $(-p_delta)")
        for (i,load) in load_active
            load["pd"] -= p_delta
        end
    end

    time_contingencies = time() - time_contingencies_start
    info(_LOGGER, "contingency eval time: $(time_contingencies)")

    return (gen_cuts=gen_cuts, branch_cuts=branch_cuts)
end




""
function check_c1_contingencies_branch_power_bpv_remote(cont_range, output_dir, cut_limit=1, solution_file="solution1.txt")
    if length(c1_network_global) <= 0 || length(c1_contingency_order_global) <= 0
        error(_LOGGER, "check_contingencies_branch_flow_remote called before load_c1_network_global")
    end

    sol = read_c1_solution1(c1_network_global, output_dir=output_dir, state_file=solution_file)
    _PM.update_data!(c1_network_global, sol)

    active_cuts = read_c1_active_flow_cuts(output_dir=output_dir)
    gen_flow_cuts = []
    branch_flow_cuts = []
    for cut in active_cuts
        if cut.cont_type == "gen"
            push!(gen_flow_cuts, cut)
        elseif cut.cont_type == "branch"
            push!(branch_flow_cuts, cut)
        else
            warn(_LOGGER, "unknown contingency type in cut $(cut)")
        end
    end

    network = copy(c1_network_global)
    contingencies = c1_contingency_order_global[cont_range]
    network["gen_contingencies"] = [c for c in contingencies if c.type == "gen"]
    network["branch_contingencies"] = [c for c in contingencies if c.type == "branch"]

    cuts = check_c1_contingencies_branch_power_bpv(network, total_cut_limit=cut_limit, gen_flow_cuts=gen_flow_cuts, branch_flow_cuts=branch_flow_cuts)

    return cuts
end


"""
A variant of `check_contingencies_branch_flow_remote`, which ignores the
participation factor based generator response model from the ARPA-e GOC
Challenge 1 specification and instead injects active power equally at all buses
in the network.
"""
function check_c1_contingencies_branch_power_bpv(network;
        gen_flow_cut_limit=10, branch_flow_cut_limit=10, total_cut_limit=typemax(Int64),
        gen_eval_limit=typemax(Int64), branch_eval_limit=typemax(Int64), sm_threshold=0.01,
        gen_flow_cuts=[], branch_flow_cuts=[]
        )

    if _IM.ismultinetwork(network)
        error(_LOGGER, "the branch flow cut generator can only be used on single networks")
    end
    time_contingencies_start = time()

    gen_cuts_active = Dict()
    for gen_cut in gen_flow_cuts
        if !haskey(gen_cuts_active, gen_cut.cont_label)
            gen_cuts_active[gen_cut.cont_label] = Set{Int}()
        end
        push!(gen_cuts_active[gen_cut.cont_label], gen_cut.branch_id)
    end

    branch_cuts_active = Dict()
    for branch_cut in branch_flow_cuts
        if !haskey(branch_cuts_active, branch_cut.cont_label)
            branch_cuts_active[branch_cut.cont_label] = Set{Int}()
        end
        push!(branch_cuts_active[branch_cut.cont_label], branch_cut.branch_id)
    end


    network_lal = deepcopy(network) #lal -> losses as loads
    ref_bus_id = _PM.reference_bus(network_lal)["index"]
    network_lal["delta"] = 0.0

    max_gen_id = maximum(parse(Int, i) for (i,gen) in network_lal["gen"])
    gen_id = prod(10 for i in 1:(Int(ceil(log10(max_gen_id))) + 1))

    gen_phantom = Dict{String,Any}()
    for (i,bus) in network_lal["bus"]
        if bus["bus_type"] != 4
            gen = Dict(
                "index" => gen_id,
                "gen_status" => 1,
                "gen_bus" => bus["index"],
                "pmin" => -Inf,
                "pmax" => Inf,
                "qmin" => -Inf,
                "qmax" => Inf,
                "pg" => 0.0,
                "qg" => 0.0
            )

            gen_phantom["$(gen_id)"] = gen
            network_lal["gen"]["$(gen_id)"] = gen
            gen_id += 1
        end
    end
    #println(keys(gen_phantom))

    load_active = Dict(i => load for (i,load) in network_lal["load"] if load["status"] != 0)

    pd_total = sum(load["pd"] for (i,load) in load_active)
    p_losses = sum(gen["pg"] for (i,gen) in network_lal["gen"] if gen["gen_status"] != 0) - pd_total
    p_delta = 0.0

    if p_losses > C1_PG_LOSS_TOL
        load_count = length(load_active)
        p_delta = p_losses/load_count
        for (i,load) in load_active
            load["pd"] += p_delta
        end
        warn(_LOGGER, "active power losses found $(p_losses) increasing loads by $(p_delta)")
    end


    gen_contingencies = calc_c1_gen_contingency_subset(network_lal, gen_eval_limit=gen_eval_limit)
    branch_contingencies = calc_c1_branch_contingency_subset(network_lal, branch_eval_limit=branch_eval_limit)

    gen_cuts = []
    for (i,cont) in enumerate(gen_contingencies)
        if length(gen_cuts) >= gen_flow_cut_limit
            info(_LOGGER, "hit gen flow cut limit $(gen_flow_cut_limit)")
            break
        end
        if length(gen_cuts) >= total_cut_limit
            info(_LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end
        #info(_LOGGER, "working on ($(i)/$(gen_eval_limit)/$(gen_cont_total)): $(cont.label)")

        cont_gen = network_lal["gen"]["$(cont.idx)"]
        pg_lost = cont_gen["pg"]

        cont_gen["gen_status"] = 0
        cont_gen["pg"] = 0.0

        pg_offset = pg_lost/length(gen_phantom)
        for (i,gen) in gen_phantom
            gen["pg"] = pg_offset
        end

        try
            solution = _PM.compute_dc_pf(network_lal)["solution"]
            _PM.update_data!(network_lal, solution)
        catch exception
            warn(_LOGGER, "linear solve failed on $(cont.label)")
            continue
        end

        flow = _PM.calc_branch_flow_dc(network_lal)
        _PM.update_data!(network_lal, flow)


        vio = calc_c1_violations(network_lal, network_lal)

        #info(_LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vios = branch_c1_violations_sorted(network_lal, network_lal)
            branch_vio = branch_vios[1]

            if !haskey(gen_cuts_active, cont.label) || !(branch_vio.branch_id in gen_cuts_active[cont.label])
                info(_LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                am = _PM.calc_susceptance_matrix(network_lal)
                branch = network_lal["branch"]["$(branch_vio.branch_id)"]

                bus_injection = calc_c1_branch_ptdf_single(am, ref_bus_id, branch)
                cut = (gen_id=cont.idx, cont_label=cont.label, branch_id=branch_vio.branch_id, rating_level=1.0, bus_injection=bus_injection)
                push!(gen_cuts, cut)
            else
                warn(_LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end
        end

        cont_gen["gen_status"] = 1
        cont_gen["pg"] = pg_lost
    end

    for (i,gen) in gen_phantom
        gen["pg"] = 0.0
    end

    # work around for julia compiler bug
    #num_buses = length(network_ref[:bus])
    #b_cont = zeros(Float64, num_buses, num_buses)

    branch_cuts = []
    for (i,cont) in enumerate(branch_contingencies)
        if length(branch_cuts) >= branch_flow_cut_limit
            info(_LOGGER, "hit branch flow cut limit $(branch_flow_cut_limit)")
            break
        end
        if length(gen_cuts) + length(branch_cuts) >= total_cut_limit
            info(_LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end

        #info(_LOGGER, "working on ($(i)/$(branch_eval_limit)/$(branch_cont_total)): $(cont.label)")

        cont_branch = network_lal["branch"]["$(cont.idx)"]
        cont_branch["br_status"] = 0

        try
            solution = _PM.compute_dc_pf(network_lal)["solution"]
            _PM.update_data!(network_lal, solution)
        catch exception
            warn(_LOGGER, "linear solve failed on $(cont.label)")
            continue
        end

        flow = _PM.calc_branch_flow_dc(network_lal)
        _PM.update_data!(network_lal, flow)

        vio = calc_c1_violations(network_lal, network_lal)

        #info(_LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vio = branch_c1_violations_sorted(network_lal, network_lal)[1]
            if !haskey(branch_cuts_active, cont.label) || !(branch_vio.branch_id in branch_cuts_active[cont.label])
                info(_LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                am = _PM.calc_susceptance_matrix(network_lal)
                branch = network_lal["branch"]["$(branch_vio.branch_id)"]

                bus_injection = calc_c1_branch_ptdf_single(am, ref_bus_id, branch)
                cut = (cont_label=cont.label, branch_id=branch_vio.branch_id, rating_level=1.0, bus_injection=bus_injection)
                push!(branch_cuts, cut)
            else
                warn(_LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end
        end

        cont_branch["br_status"] = 1
    end

    if p_delta != 0.0
        warn(_LOGGER, "re-adjusting loads by $(-p_delta)")
        for (i,load) in load_active
            load["pd"] -= p_delta
        end
    end

    time_contingencies = time() - time_contingencies_start
    info(_LOGGER, "contingency eval time: $(time_contingencies)")

    return (gen_cuts=gen_cuts, branch_cuts=branch_cuts)
end


