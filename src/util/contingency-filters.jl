
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

    cuts = check_contingencies_branch_flow_ratec_nd_first_lazy(network, total_cut_limit=cut_limit, gen_flow_cuts=gen_flow_cuts, branch_flow_cuts=branch_flow_cuts)

    return cuts
end

"given a network, checks the operating point against the contingencies to look for violations"
function check_contingencies_branch_flow_ratec_nd_first_lazy(network;
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


