
function compute_inverse_size(network)
    b0 = Base.gc_bytes()
    network_ref = PowerModels.build_ref(network)
    network_ref = network_ref[:nw][0]
    bus_idx2id, bus_id2idx = build_indices(network_ref[:bus])
    b_base = compute_B_sparse(network_ref, bus_idx2id, bus_id2idx)
    b_inv = compute_B_inverse(network_ref, bus_idx2id, bus_id2idx)
    return Base.gc_bytes() - b0
end


"maps dict ids into the range 1-to-n"
function build_indices(d::Dict)
    original_indices = sort(collect(keys(d)));
    num_indices = length(original_indices)
    reverse_indices = Dict(zip(original_indices,1:num_indices))
    return original_indices, reverse_indices
end


"computes the bus admittance matrix B, with indices given by `bus_id2idx::Dict`"
function compute_B_sparse(ref::Dict{Symbol,<:Any}, bus_idx2id::Vector{Int}, bus_id2idx::Dict{Int,Int}; inactive_branches::Set{Int}=Set(Int[]))
    num_buses = length(ref[:bus])
    I = Int64[]
    J = Int64[]
    V = Float64[]

    for (i,branch) in ref[:branch]
        if !(branch["index"] in inactive_branches)
            f_bus = bus_id2idx[ref[:bus][branch["f_bus"]]["index"]]
            t_bus = bus_id2idx[ref[:bus][branch["t_bus"]]["index"]]
            b_val = -branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2)
            push!(I, f_bus); push!(J, t_bus); push!(V,  b_val)
            push!(I, t_bus); push!(J, f_bus); push!(V,  b_val)
            push!(I, f_bus); push!(J, f_bus); push!(V, -b_val)
            push!(I, t_bus); push!(J, t_bus); push!(V, -b_val)
        end
    end

    return sparse(I,J,V)
end


"computes the bus admittance matrix B inverse, with indices given by `bus_id2idx::Dict`"
function compute_B_inverse(ref::Dict{Symbol,<:Any}, bus_idx2id::Vector{Int}, bus_id2idx::Dict{Int,Int}; inactive_branches::Set{Int}=Set(Int[]))
    num_buses = length(ref[:bus])
    B = compute_B(ref, bus_idx2id, bus_id2idx; inactive_branches=inactive_branches)

    ref_bus_id, ref_bus = first(ref[:ref_buses])
    r = bus_id2idx[ref_bus["index"]]
    nonref_indices = Int64[b for b in 1:num_buses if b != r]
    inv_B = zeros(Float64, num_buses, num_buses)
    inv_B[nonref_indices, nonref_indices] = inv(B[nonref_indices, nonref_indices])

    # matrix is dence, does not seem to help performance
    #inv_B = sparse(inv_B)
    #println(inv_B)

    return inv_B
end


"computes the bus admittance matrix B, with indices given by `bus_id2idx::Dict`"
function compute_B(ref::Dict{Symbol,<:Any}, bus_idx2id::Vector{Int}, bus_id2idx::Dict{Int,Int}; inactive_branches::Set{Int}=Set(Int[]))
    num_buses = length(ref[:bus])
    B = zeros(Float64, num_buses, num_buses)

    for (i,branch) in ref[:branch]
        if !(branch["index"] in inactive_branches)
            f_bus = bus_id2idx[ref[:bus][branch["f_bus"]]["index"]]
            t_bus = bus_id2idx[ref[:bus][branch["t_bus"]]["index"]]
            b_val = -branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2)
            B[f_bus, t_bus] += b_val
            B[t_bus, f_bus] += b_val
            B[f_bus, f_bus] += -b_val
            B[t_bus, t_bus] += -b_val
        end
    end

    return B
end


"compute bus injections"
function compute_bus_injections(ref::Dict{Symbol,<:Any}, bus_idx2id::Vector{Int})
    @assert length(ref[:bus]) == length(bus_idx2id)
    num_buses = length(bus_idx2id)
    bus_inj = [0.0 for i in 1:num_buses]
    for i in 1:num_buses
        for g in ref[:bus_gens][bus_idx2id[i]]
            bus_inj[i] += ref[:gen][g]["pg"]
        end
        for l in ref[:bus_loads][bus_idx2id[i]]
            bus_inj[i] -= ref[:load][l]["pd"]
        end
        for s in ref[:bus_shunts][bus_idx2id[i]]
            bus_inj[i] -= ref[:shunt][s]["gs"]*1.0^2
        end
    end
    return bus_inj
end

"compute bus injections"
function compute_bus_injections(ref::Dict{Symbol,<:Any}, gen_setpoint::Dict{String,<:Any}, bus_idx2id::Vector{Int})
    @assert length(ref[:bus]) == length(bus_idx2id)
    num_buses = length(bus_idx2id)
    bus_inj = [0.0 for i in 1:num_buses]
    for i in 1:num_buses
        for g in ref[:bus_gens][bus_idx2id[i]]
            bus_inj[i] += gen_setpoint["$(g)"]["pg"]
        end
        for l in ref[:bus_loads][bus_idx2id[i]]
            bus_inj[i] -= ref[:load][l]["pd"]
        end
        for s in ref[:bus_shunts][bus_idx2id[i]]
            bus_inj[i] -= ref[:shunt][s]["gs"]*1.0^2
        end
    end
    return bus_inj
end

"compute bus injections"
function compute_bus_injections(ref::Dict{Symbol,<:Any}, gen_setpoint::Dict{String,<:Any}, load_setpoint::Dict{String,<:Any}, bus_idx2id::Vector{Int}; bus_inj_offset::Float64=0.0)
    @assert length(ref[:bus]) == length(bus_idx2id)
    num_buses = length(bus_idx2id)
    bus_inj = [bus_inj_offset for i in 1:num_buses]
    for i in 1:num_buses
        for g in ref[:bus_gens][bus_idx2id[i]]
            bus_inj[i] += gen_setpoint["$(g)"]["pg"]
        end
        for l in ref[:bus_loads][bus_idx2id[i]]
            bus_inj[i] -= load_setpoint["$(l)"]["pd"]
        end
        for s in ref[:bus_shunts][bus_idx2id[i]]
            bus_inj[i] -= ref[:shunt][s]["gs"]*1.0^2
        end
    end
    return bus_inj
end


function solve_theta(B, bus_injection::Vector{Float64}, ref_bus::Int64)
    theta = B \ bus_injection
    theta = theta .- theta[ref_bus]
    return theta
end

compute_flows(ref::Dict{Symbol,<:Any}, theta, bus_id2idx::Dict{Int,Int}) = Dict(l => (-branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))*(theta[bus_id2idx[branch["f_bus"]]] - theta[bus_id2idx[branch["t_bus"]]]) for (l,branch) in ref[:branch])


function compute_branch_ptdf_single(ref::Dict{Symbol,<:Any}, B_inverse::Matrix{Float64}, bus_id2idx::Dict{Int,Int}, branch_id::Int)
    branch_ptdf = Dict{Int,Any}()
    branch = ref[:branch][branch_id]
    row_fr = bus_id2idx[branch["f_bus"]]
    row_to = bus_id2idx[branch["t_bus"]]
    b = (-branch["br_x"]/(branch["br_x"]^2+branch["br_r"]^2))

    branch_ptdf[branch_id] = -b*(B_inverse[row_fr,:] - B_inverse[row_to,:])
    return branch_ptdf
end

function built_flow_cut(cont_label, branch_id, rating_level, branch_ptdf::Dict{Int,<:Any}, bus_idx2id::Vector{Int})
    bus_injection = Dict(bus_idx2id[idx] => value for (idx,value) in enumerate(branch_ptdf[branch_id]))
    cut = (cont_label=cont_label, branch_id=branch_id, rating_level=rating_level, bus_injection=bus_injection)
end



# a global network variable used for iterative computations
network_global = Dict{String,Any}()

# a global contingency list used for iterative computations
contingency_order_global = []

""
function load_network_global(con_file, inl_file, raw_file, rop_file, scenario_id)
    info(LOGGER, "skipping goc and power models data warnings")
    pm_logger_level = getlevel(getlogger(PowerModels))
    goc_logger_level = getlevel(LOGGER)

    setlevel!(getlogger(PowerModels), "error")
    setlevel!(LOGGER, "error")

    goc_data = parse_goc_files(con_file, inl_file, raw_file, rop_file, scenario_id=scenario_id)
    global network_global = build_pm_model(goc_data)
    global contingency_order_global = contingency_order(network_global)

    setlevel!(getlogger(PowerModels), pm_logger_level)
    setlevel!(LOGGER, goc_logger_level)

    return 0
end


"build a static ordering of all contigencies"
function contingency_order(network)
    gen_cont_order = sort(network["gen_contingencies"], by=(x) -> x.label)
    branch_cont_order = sort(network["branch_contingencies"], by=(x) -> x.label)

    gen_cont_total = length(gen_cont_order)
    branch_cont_total = length(branch_cont_order)

    gen_rate = 1.0
    branch_rate = 1.0
    steps = 1

    if gen_cont_total == 0 && branch_cont_total == 0
        # defaults are good
    elseif gen_cont_total == 0 && branch_cont_total != 0
        steps = branch_cont_total
    elseif gen_cont_total != 0 < branch_cont_total == 0
        steps = gen_cont_total
    elseif gen_cont_total == branch_cont_total
        steps = branch_cont_total
    elseif gen_cont_total < branch_cont_total
        gen_rate = 1.0
        branch_rate = branch_cont_total/gen_cont_total
        steps = gen_cont_total
    elseif gen_cont_total > branch_cont_total
        gen_rate = gen_cont_total/branch_cont_total
        branch_rate = 1.0 
        steps = branch_cont_total
    end


    #println(gen_cont_total)
    #println(branch_cont_total)
    #println(steps)

    #println(gen_rate)
    #println(branch_rate)
    #println("")

    cont_order = []
    gen_cont_start = 1
    branch_cont_start = 1
    for s in 1:steps
        gen_cont_end = min(gen_cont_total, trunc(Int,ceil(s*gen_rate)))
        #println(gen_cont_start:gen_cont_end)
        for j in gen_cont_start:gen_cont_end
            push!(cont_order, gen_cont_order[j])
        end
        gen_cont_start = gen_cont_end+1

        branch_cont_end = min(branch_cont_total, trunc(Int,ceil(s*branch_rate)))
        #println("$(s) - $(branch_cont_start:branch_cont_end)")
        for j in branch_cont_start:branch_cont_end
            push!(cont_order, branch_cont_order[j])
        end
        branch_cont_start = branch_cont_end+1
    end

    #=
    for s in 1:steps
        gen_cont_start = trunc(Int, ceil(1+(s-1)*gen_rate))
        gen_cont_end = min(gen_cont_total, trunc(Int,ceil(s*gen_rate)))
        #println(gen_cont_start:gen_cont_end)
        for j in gen_cont_start:gen_cont_end
            push!(cont_order, gen_cont_order[j])
        end

        branch_cont_start = trunc(Int, ceil(1+(s-1)*branch_rate))
        branch_cont_end = min(branch_cont_total, trunc(Int,ceil(s*branch_rate)))
        println("$(s) - $(branch_cont_start:branch_cont_end)")
        for j in branch_cont_start:branch_cont_end
            push!(cont_order, branch_cont_order[j])
        end
    end
    =#

    #println(length(cont_order))
    #println(gen_cont_total + branch_cont_total)

    @assert(length(cont_order) == gen_cont_total + branch_cont_total)

    return cont_order
end


function write_contingencies(network; output_dir="", file_name="contingencies.txt")
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

function write_active_flow_cuts(network; output_dir="", file_name="active_flow_cuts.txt")
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

function read_active_flow_cuts(;output_dir="", file_name="active_flow_cuts.txt")
    if length(output_dir) > 0
        file_path = joinpath(output_dir, file_name)
    else
        file_path = file_name
    end

    if isfile(file_path)
        info(LOGGER, "loading flow cuts file: $(file_path)")
        return parse_flow_cuts_file(file_path)
    else
        info(LOGGER, "flow cuts file not found: $(file_path)")
        return []
    end
end

function parse_flow_cuts_file(file::String)
    open(file) do io
        return parse_flow_cuts_file(io)
    end
end

function parse_flow_cuts_file(io::IO)
    cuts_list = []

    for line in readlines(io)
        if length(strip(line)) == 0
            warn(LOGGER, "skipping blank line in cuts file")
            continue
        end
        line_parts = split(line, ",")
        if length(line_parts) != 4
            warn(LOGGER, "skipping ill formated line\n   $(line)")
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



""
function check_contingencies_branch_flow_remote_nd_first_lazy(cont_range, output_dir, cut_limit=1, solution_file="solution1.txt")
    if length(network_global) <= 0 || length(contingency_order_global) <= 0
        error(LOGGER, "check_contingencies_branch_flow_remote called before load_network_global")
    end

    sol = read_solution1(network_global, output_dir=output_dir, state_file=solution_file)
    PowerModels.update_data!(network_global, sol)

    active_cuts = read_active_flow_cuts(output_dir=output_dir)
    gen_flow_cuts = []
    branch_flow_cuts = []
    for cut in active_cuts
        if cut.cont_type == "gen"
            push!(gen_flow_cuts, cut)
        elseif cut.cont_type == "branch"
            push!(branch_flow_cuts, cut)
        else
            warn(LOGGER, "unknown contingency type in cut $(cut)")
        end
    end

    network = copy(network_global)
    contingencies = contingency_order_global[cont_range]
    network["gen_contingencies"] = [c for c in contingencies if c.type == "gen"]
    network["branch_contingencies"] = [c for c in contingencies if c.type == "branch"]

    cuts = check_contingencies_branch_power_bpv(network, total_cut_limit=cut_limit, gen_flow_cuts=gen_flow_cuts, branch_flow_cuts=branch_flow_cuts)

    return cuts
end


"""
A variant of check_contingencies_branch_flow_remote, which ignores the
participation factor based generator response model from the ARPA-e GOC
Challenge 1 specification and instead injects active power equally at all buses
in the network.
"""
function check_contingencies_branch_power_bpv(network;
        gen_flow_cut_limit=10, branch_flow_cut_limit=10, total_cut_limit=typemax(Int64),
        gen_eval_limit=typemax(Int64), branch_eval_limit=typemax(Int64), sm_threshold=0.01,
        gen_flow_cuts=[], branch_flow_cuts=[]
        )

    if InfrastructureModels.ismultinetwork(network)
        error(LOGGER, "the branch flow cut generator can only be used on single networks")
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


    network_ref = PowerModels.build_ref(network)
    network_ref = network_ref[:nw][0]


    solution_base = extract_solution(network; branch_flow=true)

    pd_total = sum(load["pd"] for (i,load) in network_ref[:load])
    p_losses = sum(gen["pg"] for (i,gen) in network_ref[:gen]) - pd_total
    load_setpoint = network["load"]
    p_delta = 0.0

    if p_losses > pg_loss_tol
        load_count = length(load_setpoint)
        p_delta = p_losses/load_count
        for (i,load) in load_setpoint
            load["pd"] = load["pd"] + p_delta
        end
        warn(LOGGER, "active power losses found $(p_losses) increasing loads by $(p_delta)")
    end


    bus_idx2id, bus_id2idx = build_indices(network_ref[:bus])
    b_base = compute_B_sparse(network_ref, bus_idx2id, bus_id2idx)
    b_inv_base = compute_B_inverse(network_ref, bus_idx2id, bus_id2idx)
    #branch_ptdf_base = compute_branch_ptdf(network_ref, b_inv_base, bus_id2idx)
    bus_inj_base = compute_bus_injections(network_ref, solution_base["gen"], load_setpoint, bus_idx2id)

    @assert length(network_ref[:ref_buses]) == 1
    ref_bus_id, ref_bus = first(network_ref[:ref_buses])
    ref_bus_idx = bus_id2idx[ref_bus["index"]]

    gen_cont_total = length(network["gen_contingencies"])
    branch_cont_total = length(network["branch_contingencies"])

    gen_eval_limit = min(gen_eval_limit, gen_cont_total)
    branch_eval_limit = min(branch_eval_limit, branch_cont_total)

    gen_cap = Dict(gen["index"] => sqrt(max(abs(gen["pmin"]), abs(gen["pmax"]))^2 + max(abs(gen["qmin"]), abs(gen["qmax"]))^2) for (i,gen) in network["gen"])
    network["gen_contingencies"] = sort(network["gen_contingencies"], rev=true, by=x -> gen_cap[x.idx])
    gen_contingencies = network["gen_contingencies"][1:gen_eval_limit]

    line_imp_mag = Dict(branch["index"] => branch["rate_a"]*sqrt(branch["br_r"]^2 + branch["br_x"]^2) for (i,branch) in network["branch"])
    network["branch_contingencies"] = sort(network["branch_contingencies"], rev=true, by=x -> line_imp_mag[x.idx])
    branch_contingencies = network["branch_contingencies"][1:branch_eval_limit]

    gen_cuts = []
    for (i,cont) in enumerate(gen_contingencies)
        if length(gen_cuts) >= gen_flow_cut_limit
            info(LOGGER, "hit gen flow cut limit $(gen_flow_cut_limit)")
            break
        end
        if length(gen_cuts) >= total_cut_limit
            info(LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end
        #info(LOGGER, "working on ($(i)/$(gen_eval_limit)/$(gen_cont_total)): $(cont.label)")

        sol_tmp = deepcopy(solution_base)
        cont_gen = sol_tmp["gen"]["$(cont.idx)"]
        pg_lost = cont_gen["pg"]
        qg_lost = cont_gen["qg"]
        cont_gen["pg"] = 0.0
        cont_gen["qg"] = 0.0

        gen = network_ref[:gen][cont.idx]
        gen_bus = network_ref[:bus][gen["gen_bus"]]

        sol_tmp["delta"] = 0.0

        bus_inj_cont = compute_bus_injections(network_ref, sol_tmp["gen"], load_setpoint, bus_idx2id, bus_inj_offset=pg_lost/length(bus_idx2id))

        va_vector = []
        try
            va_vector = solve_theta(b_base, bus_inj_cont, ref_bus_idx)
        catch exception
            warn(LOGGER, "linear solve failed on $(cont.label)")
            continue
        end
        p_flows = compute_flows(network_ref, va_vector, bus_id2idx)

        for (l,flow) in p_flows
            sol_branch = sol_tmp["branch"]["$(l)"]
            sol_branch["pf"] =  flow
            sol_branch["pt"] = -flow
            sol_branch["qf"] = 0.0
            sol_branch["qt"] = 0.0
        end

        #PowerModels.print_summary(sol_tmp)

        network["gen"]["$(cont.idx)"]["gen_status"] = 0
        vio = compute_violations_ratec(network, sol_tmp)
        network["gen"]["$(cont.idx)"]["gen_status"] = 1

        #info(LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vio = branch_violations_sorted_ratec(network, sol_tmp)[1]
            if !haskey(gen_cuts_active, cont.label) || !(branch_vio.branch_id in gen_cuts_active[cont.label])
                info(LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                branch_ptdf_cont = compute_branch_ptdf_single(network_ref, b_inv_base, bus_id2idx, branch_vio.branch_id)
                cut = built_flow_cut(cont.label, branch_vio.branch_id, 1.0, branch_ptdf_cont, bus_idx2id)
                cut = (gen_id=cont.idx, cont_label=cut.cont_label, branch_id=cut.branch_id, rating_level=cut.rating_level, bus_injection=cut.bus_injection)
                push!(gen_cuts, cut)
                break
            else
                warn(LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end
        end
    end

    # work around for julia compiler bug
    #num_buses = length(network_ref[:bus])
    #b_cont = zeros(Float64, num_buses, num_buses)

    branch_cuts = []
    for (i,cont) in enumerate(branch_contingencies)
        if length(branch_cuts) >= branch_flow_cut_limit
            info(LOGGER, "hit branch flow cut limit $(branch_flow_cut_limit)")
            break
        end
        if length(gen_cuts) + length(branch_cuts) >= total_cut_limit
            info(LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end

        #info(LOGGER, "working on ($(i)/$(branch_eval_limit)/$(branch_cont_total)): $(cont.label)")
        sol_tmp = deepcopy(solution_base)
        if haskey(sol_tmp, "branch")
            cont_branch = sol_tmp["branch"]["$(cont.idx)"]
            cont_branch["pf"] = 0.0
            cont_branch["pt"] = 0.0
            cont_branch["qf"] = 0.0
            cont_branch["qt"] = 0.0
        end

        sol_tmp["delta"] = 0.0

        b_cont = compute_B_sparse(network_ref, bus_idx2id, bus_id2idx, inactive_branches=Set([cont.idx]))
        # work around for julia compiler bug
        #update_B!(b_cont, network_ref, bus_idx2id, bus_id2idx, inactive_branches=Set([cont.idx]))

        va_vector = []
        try
            va_vector = solve_theta(b_cont, bus_inj_base, ref_bus_idx)
        catch exception
            warn(LOGGER, "linear solve failed on $(cont.label)")
            continue
        end
        p_flows = compute_flows(network_ref, va_vector, bus_id2idx)

        for (l,flow) in p_flows
            sol_branch = sol_tmp["branch"]["$(l)"]
            if l != cont.idx
                sol_branch["pf"] =  flow
                sol_branch["pt"] = -flow
            else
                sol_branch["pf"] = 0.0
                sol_branch["pt"] = 0.0
            end

            sol_branch["qf"] = 0.0
            sol_branch["qt"] = 0.0
        end

        #PowerModels.print_summary(sol_tmp)

        network["branch"]["$(cont.idx)"]["br_status"] = 0
        vio = compute_violations_ratec(network, sol_tmp)
        network["branch"]["$(cont.idx)"]["br_status"] = 1

        #info(LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vio = branch_violations_sorted_ratec(network, sol_tmp)[1]

            if !haskey(branch_cuts_active, cont.label) || !(branch_vio.branch_id in branch_cuts_active[cont.label])
                info(LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                b_inv_cont = compute_B_inverse(network_ref, bus_idx2id, bus_id2idx, inactive_branches=Set([cont.idx]))
                branch_ptdf_cont = compute_branch_ptdf_single(network_ref, b_inv_cont, bus_id2idx, branch_vio.branch_id)

                cut = built_flow_cut(cont.label, branch_vio.branch_id, 1.0, branch_ptdf_cont, bus_idx2id)
                push!(branch_cuts, cut)
                break
            else
                warn(LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end
        end
    end

    if p_delta != 0.0
        warn(LOGGER, "re-adjusting loads by $(-p_delta)")
        for (i,load) in load_setpoint
            load["pd"] = load["pd"] - p_delta
        end
    end

    time_contingencies = time() - time_contingencies_start
    info(LOGGER, "cont eval time: $(time_contingencies)")

    return (gen_cuts=gen_cuts, branch_cuts=branch_cuts)
end








""
function check_contingencies_branch_flow_remote(cont_range, output_dir, cut_limit=1, solution_file="solution1.txt")
    if length(network_global) <= 0 || length(contingency_order_global) <= 0
        error(LOGGER, "check_contingencies_branch_flow_remote called before load_network_global")
    end

    sol = read_solution1(network_global, output_dir=output_dir, state_file=solution_file)
    PowerModels.update_data!(network_global, sol)

    active_cuts = read_active_flow_cuts(output_dir=output_dir)
    gen_flow_cuts = []
    branch_flow_cuts = []
    for cut in active_cuts
        if cut.cont_type == "gen"
            push!(gen_flow_cuts, cut)
        elseif cut.cont_type == "branch"
            push!(branch_flow_cuts, cut)
        else
            warn(LOGGER, "unknown contingency type in cut $(cut)")
        end
    end

    network = copy(network_global)
    contingencies = contingency_order_global[cont_range]
    network["gen_contingencies"] = [c for c in contingencies if c.type == "gen"]
    network["branch_contingencies"] = [c for c in contingencies if c.type == "branch"]

    cuts = check_contingencies_branch_power(network, total_cut_limit=cut_limit, gen_flow_cuts=gen_flow_cuts, branch_flow_cuts=branch_flow_cuts)

    return cuts
end



"""
Checks a given operating point against the contingencies to look for branch
flow violations.  The DC Power Flow approximation is used for flow simulation.
If a violation is found, computes a PTDF cut based on bus injections.  Uses the
participation factor based generator response model from the ARPA-e GOC
Challenge 1 specification.
"""
function check_contingencies_branch_power(network;
        gen_flow_cut_limit=10, branch_flow_cut_limit=10, total_cut_limit=typemax(Int64),
        gen_eval_limit=typemax(Int64), branch_eval_limit=typemax(Int64), sm_threshold=0.01,
        gen_flow_cuts=[], branch_flow_cuts=[]
        )

    if InfrastructureModels.ismultinetwork(network)
        error(LOGGER, "the branch flow cut generator can only be used on single networks")
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


    network_ref = PowerModels.build_ref(network)
    network_ref = network_ref[:nw][0]


    solution_base = extract_solution(network; branch_flow=true)

    pd_total = sum(load["pd"] for (i,load) in network_ref[:load])
    p_losses = sum(gen["pg"] for (i,gen) in network_ref[:gen]) - pd_total
    load_setpoint = network["load"]
    p_delta = 0.0

    if p_losses > pg_loss_tol
        load_count = length(load_setpoint)
        p_delta = p_losses/load_count
        for (i,load) in load_setpoint
            load["pd"] = load["pd"] + p_delta
        end
        warn(LOGGER, "active power losses found $(p_losses) increasing loads by $(p_delta)")
    end


    bus_idx2id, bus_id2idx = build_indices(network_ref[:bus])
    b_base = compute_B_sparse(network_ref, bus_idx2id, bus_id2idx)
    b_inv_base = compute_B_inverse(network_ref, bus_idx2id, bus_id2idx)
    #branch_ptdf_base = compute_branch_ptdf(network_ref, b_inv_base, bus_id2idx)
    bus_inj_base = compute_bus_injections(network_ref, solution_base["gen"], load_setpoint, bus_idx2id)

    @assert length(network_ref[:ref_buses]) == 1
    ref_bus_id, ref_bus = first(network_ref[:ref_buses])
    ref_bus_idx = bus_id2idx[ref_bus["index"]]

    gen_cont_total = length(network["gen_contingencies"])
    branch_cont_total = length(network["branch_contingencies"])

    gen_eval_limit = min(gen_eval_limit, gen_cont_total)
    branch_eval_limit = min(branch_eval_limit, branch_cont_total)

    gen_cap = Dict(gen["index"] => sqrt(max(abs(gen["pmin"]), abs(gen["pmax"]))^2 + max(abs(gen["qmin"]), abs(gen["qmax"]))^2) for (i,gen) in network["gen"])
    network["gen_contingencies"] = sort(network["gen_contingencies"], rev=true, by=x -> gen_cap[x.idx])
    gen_contingencies = network["gen_contingencies"][1:gen_eval_limit]

    line_imp_mag = Dict(branch["index"] => branch["rate_a"]*sqrt(branch["br_r"]^2 + branch["br_x"]^2) for (i,branch) in network["branch"])
    network["branch_contingencies"] = sort(network["branch_contingencies"], rev=true, by=x -> line_imp_mag[x.idx])
    branch_contingencies = network["branch_contingencies"][1:branch_eval_limit]

    gen_cuts = []
    for (i,cont) in enumerate(gen_contingencies)
        if length(gen_cuts) >= gen_flow_cut_limit
            info(LOGGER, "hit gen flow cut limit $(gen_flow_cut_limit)")
            break
        end
        if length(gen_cuts) >= total_cut_limit
            info(LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end
        #info(LOGGER, "working on ($(i)/$(gen_eval_limit)/$(gen_cont_total)): $(cont.label)")

        sol_tmp = deepcopy(solution_base)
        cont_gen = sol_tmp["gen"]["$(cont.idx)"]
        pg_lost = cont_gen["pg"]
        qg_lost = cont_gen["qg"]
        cont_gen["pg"] = 0.0
        cont_gen["qg"] = 0.0

        gen = network_ref[:gen][cont.idx]
        gen_bus = network_ref[:bus][gen["gen_bus"]]
        gen_set = network_ref[:area_gens][gen_bus["area"]]

        alpha_gens = [gen["alpha"] for (i,gen) in network_ref[:gen] if gen["index"] != cont.idx && gen["index"] in gen_set]
        if length(alpha_gens) == 0 || isapprox(sum(alpha_gens), 0.0)
            warn(LOGGER, "no available active power response in cont $(cont.label), active gens $(length(alpha_gens))")
            continue
        end
        alpha_total = sum(alpha_gens)
        delta = pg_lost/alpha_total
        sol_tmp["delta"] = delta
        #info(LOGGER, "$(pg_lost) - $(alpha_total) - $(delta)")

        for (i,gen) in network_ref[:gen]
            sol_gen = sol_tmp["gen"]["$(i)"]
            if gen["index"] != cont.idx && gen["index"] in gen_set
                sol_gen["pg"] = sol_gen["pg"] + gen["alpha"]*delta
            end
            sol_gen["qg"] = 0.0
        end

        bus_inj_cont = compute_bus_injections(network_ref, sol_tmp["gen"], load_setpoint, bus_idx2id)
        #p_flows = compute_branch_flows(branch_ptdf_base, bus_inj_cont)
        va_vector = []
        try
            va_vector = solve_theta(b_base, bus_inj_cont, ref_bus_idx)
        catch exception
            warn(LOGGER, "linear solve failed on $(cont.label)")
            continue
        end
        p_flows = compute_flows(network_ref, va_vector, bus_id2idx)


        for (l,flow) in p_flows
            sol_branch = sol_tmp["branch"]["$(l)"]
            sol_branch["pf"] =  flow
            sol_branch["pt"] = -flow
            sol_branch["qf"] = 0.0
            sol_branch["qt"] = 0.0
        end

        #PowerModels.print_summary(sol_tmp)

        network["gen"]["$(cont.idx)"]["gen_status"] = 0
        vio = compute_violations_ratec(network, sol_tmp)
        network["gen"]["$(cont.idx)"]["gen_status"] = 1

        #info(LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vios = branch_violations_sorted_ratec(network, sol_tmp)
            branch_vio = branch_vios[1]

            if !haskey(gen_cuts_active, cont.label) || !(branch_vio.branch_id in gen_cuts_active[cont.label])
                info(LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                branch_ptdf_cont = compute_branch_ptdf_single(network_ref, b_inv_base, bus_id2idx, branch_vio.branch_id)
                cut = built_flow_cut(cont.label, branch_vio.branch_id, 1.0, branch_ptdf_cont, bus_idx2id)
                cut = (gen_id=cont.idx, cont_label=cut.cont_label, branch_id=cut.branch_id, rating_level=cut.rating_level, bus_injection=cut.bus_injection)
                push!(gen_cuts, cut)
            else
                warn(LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end

        end
    end

    # work around for julia compiler bug
    #num_buses = length(network_ref[:bus])
    #b_cont = zeros(Float64, num_buses, num_buses)

    branch_cuts = []
    for (i,cont) in enumerate(branch_contingencies)
        if length(branch_cuts) >= branch_flow_cut_limit
            info(LOGGER, "hit branch flow cut limit $(branch_flow_cut_limit)")
            break
        end
        if length(gen_cuts) + length(branch_cuts) >= total_cut_limit
            info(LOGGER, "hit total cut limit $(total_cut_limit)")
            break
        end

        #info(LOGGER, "working on ($(i)/$(branch_eval_limit)/$(branch_cont_total)): $(cont.label)")
        sol_tmp = deepcopy(solution_base)
        if haskey(sol_tmp, "branch")
            cont_branch = sol_tmp["branch"]["$(cont.idx)"]
            cont_branch["pf"] = 0.0
            cont_branch["pt"] = 0.0
            cont_branch["qf"] = 0.0
            cont_branch["qt"] = 0.0
        end

        #delta = 0.0
        for (i,gen) in network_ref[:gen]
            sol_gen = sol_tmp["gen"]["$(i)"]
            sol_gen["pg"] = sol_gen["pg"] #+ gen["alpha"]*delta
        end

        b_cont = compute_B_sparse(network_ref, bus_idx2id, bus_id2idx, inactive_branches=Set([cont.idx]))
        # work around for julia compiler bug
        #update_B!(b_cont, network_ref, bus_idx2id, bus_id2idx, inactive_branches=Set([cont.idx]))

        va_vector = []
        try
            va_vector = solve_theta(b_cont, bus_inj_base, ref_bus_idx)
        catch exception
            warn(LOGGER, "linear solve failed on $(cont.label)")
            continue
        end
        p_flows = compute_flows(network_ref, va_vector, bus_id2idx)

        for (l,flow) in p_flows
            sol_branch = sol_tmp["branch"]["$(l)"]
            if l != cont.idx
                sol_branch["pf"] =  flow
                sol_branch["pt"] = -flow
            else
                sol_branch["pf"] = 0.0
                sol_branch["pt"] = 0.0
            end

            sol_branch["qf"] = 0.0
            sol_branch["qt"] = 0.0
        end

        #PowerModels.print_summary(sol_tmp)

        network["branch"]["$(cont.idx)"]["br_status"] = 0
        vio = compute_violations_ratec(network, sol_tmp)
        network["branch"]["$(cont.idx)"]["br_status"] = 1

        #info(LOGGER, "$(cont.label) violations $(vio)")
        #if vio.vm > vm_threshold || vio.pg > pg_threshold || vio.qg > qg_threshold || vio.sm > sm_threshold
        if vio.sm > sm_threshold
            branch_vio = branch_violations_sorted_ratec(network, sol_tmp)[1]
            if !haskey(branch_cuts_active, cont.label) || !(branch_vio.branch_id in branch_cuts_active[cont.label])
                info(LOGGER, "adding flow cut on cont $(cont.label) branch $(branch_vio.branch_id) due to constraint flow violations $(branch_vio.sm_vio)")

                b_inv_cont = compute_B_inverse(network_ref, bus_idx2id, bus_id2idx, inactive_branches=Set([cont.idx]))
                branch_ptdf_cont = compute_branch_ptdf_single(network_ref, b_inv_cont, bus_id2idx, branch_vio.branch_id)

                cut = built_flow_cut(cont.label, branch_vio.branch_id, 1.0, branch_ptdf_cont, bus_idx2id)
                push!(branch_cuts, cut)
            else
                warn(LOGGER, "skipping active flow cut on cont $(cont.label) branch $(branch_vio.branch_id) with constraint flow violations $(branch_vio.sm_vio)")
            end
        end
    end

    if p_delta != 0.0
        warn(LOGGER, "re-adjusting loads by $(-p_delta)")
        for (i,load) in load_setpoint
            load["pd"] = load["pd"] - p_delta
        end
    end

    time_contingencies = time() - time_contingencies_start
    info(LOGGER, "cont eval time: $(time_contingencies)")

    return (gen_cuts=gen_cuts, branch_cuts=branch_cuts)
end

