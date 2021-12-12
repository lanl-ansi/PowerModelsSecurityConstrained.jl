

"build a static ordering of all contingencies"
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
    elseif gen_cont_total != 0 && branch_cont_total == 0
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

    @assert(length(cont_order) == gen_cont_total + branch_cont_total)

    return cont_order
end


"""
transforms a contigency list into explicit multinetwork data with network 0
being the base case
"""
function build_c1_scopf_multinetwork(network::Dict{String,<:Any})
    if _IM.ismultinetwork(network)
        error(_LOGGER, "build scopf can only be used on single networks")
    end

    contingencies = length(network["gen_contingencies"]) + length(network["branch_contingencies"])

    info(_LOGGER, "building scopf multi-network with $(contingencies+1) networks")

    if contingencies > 0
        mn_data = _PM.replicate(network, contingencies)
        base_network = mn_data["nw"]["0"] = deepcopy(mn_data["nw"]["1"])

        for (n, network) in mn_data["nw"]
            if n == "0"
                continue
            end

            for (i,bus) in network["bus"]
                if haskey(bus, "evhi")
                    bus["vmax"] = bus["evhi"]
                end
                if haskey(bus, "evlo")
                    bus["vmin"] = bus["evlo"]
                end
            end

            for (i,branch) in network["branch"]
                if haskey(branch, "rate_c")
                    branch["rate_a"] = branch["rate_c"]
                end
            end
        end

        network_id = 1
        for cont in base_network["gen_contingencies"]
            cont_nw = mn_data["nw"]["$(network_id)"]
            cont_nw["name"] = cont.label
            cont_gen = cont_nw["gen"]["$(cont.idx)"]
            cont_gen["gen_status"] = 0

            gen_buses = Set{Int}()
            for (i,gen) in cont_nw["gen"]
                if gen["gen_status"] != 0
                    push!(gen_buses, gen["gen_bus"])
                end
            end
            cont_nw["gen_buses"] = gen_buses

            network["response_gens"] = Set()
            gen_bus = cont_nw["bus"]["$(cont_gen["gen_bus"])"]
            cont_nw["response_gens"] = cont_nw["area_gens"][gen_bus["area"]]

            network_id += 1
        end
        for cont in base_network["branch_contingencies"]
            cont_nw = mn_data["nw"]["$(network_id)"]
            cont_nw["name"] = cont.label
            cont_branch = cont_nw["branch"]["$(cont.idx)"]
            cont_branch["br_status"] = 0

            gen_buses = Set{Int}()
            for (i,gen) in cont_nw["gen"]
                if gen["gen_status"] != 0
                    push!(gen_buses, gen["gen_bus"])
                end
            end
            cont_nw["gen_buses"] = gen_buses

            fr_bus = cont_nw["bus"]["$(cont_branch["f_bus"])"]
            to_bus = cont_nw["bus"]["$(cont_branch["t_bus"])"]

            cont_nw["response_gens"] = Set()
            if haskey(cont_nw["area_gens"], fr_bus["area"])
                cont_nw["response_gens"] = cont_nw["area_gens"][fr_bus["area"]]
            end
            if haskey(network["area_gens"], to_bus["area"])
                cont_nw["response_gens"] = union(cont_nw["response_gens"], cont_nw["area_gens"][to_bus["area"]])
            end

            network_id += 1
        end

    else
        mn_data = _PM.replicate(network, 1)
        mn_data["nw"]["0"] = mn_data["nw"]["1"]
        delete!(mn_data["nw"], "1")
    end

    return mn_data
end



# note this is simialr to bus_gen_lookup in PowerModels
# core differences are taking network as an arg and filtering by gen_status
function gens_by_bus(network::Dict{String,<:Any})
    bus_gens = Dict(i => Any[] for (i,bus) in network["bus"])
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            push!(bus_gens["$(gen["gen_bus"])"], gen)
        end
    end
    return bus_gens
end




function c1_tighten_constraints!(network::Dict{String,<:Any})
    for (i,bus) in network["bus"]
        if isapprox(bus["vmax"], bus["evhi"])
            bus["vmax_target"] = bus["vmax"] - 0.03
        else
            bus["vmax_target"] = bus["vmax"]
        end
        if isapprox(bus["vmin"], bus["evlo"])
            bus["vmin_target"] = bus["vmin"] + 0.03
        else
            bus["vmin_target"] = bus["vmin"]
        end
    end

    for (i,gen) in network["gen"]
        if gen["pmax"] > 0 && gen["pmax"]*0.9 > gen["pmin"]
            gen["pmax_target"] = gen["pmax"]*0.9
        else
            gen["pmax_target"] = gen["pmax"]
        end
    end

    for (i,branch) in network["branch"]
        # TODO restore this
        #if isapprox(branch["rate_a"], branch["rate_c"])
            branch["rate_a"] = branch["rate_a"]*0.80
        #end
    end
end



function deactivate_rate_a!(network::Dict{String,<:Any})
    network["active_rates"] = Int[]
    for (i,branch) in network["branch"]
        branch["rate_a_inactive"] = branch["rate_a"]
        delete!(branch, "rate_a")
    end
end

function activate_rate_a!(network::Dict{String,<:Any})
    if haskey(network, "active_rates")
        delete!(network, "active_rates")
    end

    for (i,branch) in network["branch"]
        if haskey(branch, "rate_a_inactive")
            branch["rate_a"] = branch["rate_a_inactive"]
            delete!(branch, "rate_a_inactive")
        end
    end
end

function activate_rate_a_violations!(network::Dict{String,<:Any})
    ac_flows = _PM.calc_branch_flow_ac(network)
    for (i,branch) in network["branch"]
        branch["pf_start"] = ac_flows["branch"][i]["pf"]
        branch["qf_start"] = ac_flows["branch"][i]["qf"]

        branch["pt_start"] = ac_flows["branch"][i]["pt"]
        branch["qt_start"] = ac_flows["branch"][i]["qt"]
    end

    line_flow_vio = false
    for (i,branch) in network["branch"]
        if !haskey(branch, "rate_a")
            if (ac_flows["branch"][i]["pf"]^2 + ac_flows["branch"][i]["qf"]^2 > branch["rate_a_inactive"]^2 ||
                ac_flows["branch"][i]["pt"]^2 + ac_flows["branch"][i]["qt"]^2 > branch["rate_a_inactive"]^2)
                info(_LOGGER, "add rate_a flow limit on branch $(i) $(branch["source_id"])")
                #branch["rate_a"] = branch["rate_a_inactive"] - max(abs(ac_flows["branch"][i]["qf"]), abs(ac_flows["branch"][i]["qt"]))
                branch["rate_a"] = branch["rate_a_inactive"]
                push!(network["active_rates"], branch["index"])
                line_flow_vio = true
            end
        else
            sm_fr = sqrt(ac_flows["branch"][i]["pf"]^2 + ac_flows["branch"][i]["qf"]^2)
            sm_to = sqrt(ac_flows["branch"][i]["pf"]^2 + ac_flows["branch"][i]["qf"]^2)
            vio = max(0.0, sm_fr - branch["rate_a"], sm_to - branch["rate_a"])
            if vio > 0.01
                warn(_LOGGER, "add rate_a flow limit violations $(vio) on branch $(i) $(branch["source_id"])")
            end
        end
    end

    return line_flow_vio
end


"""
assumes there is one reference bus and one connected component and adjusts voltage
angles to be centered around zero at the reference bus.
"""
function correct_voltage_angles!(network::Dict{String,<:Any})
    ref_bus = -1
    for (i,bus) in network["bus"]
        if bus["bus_type"] == 3
            @assert ref_bus == -1
            ref_bus = bus
        end
    end

    if !isapprox(ref_bus["va"], 0.0, atol=1e-8)
        warn(_LOGGER, "shifting voltage angles by $(-ref_bus["va"]) to set reference bus to 0.0")
        shift_voltage_anlges!(network, -ref_bus["va"])
    end
end


"shift networks voltage angles by a specified amount"
function shift_voltage_anlges!(network::Dict{String,<:Any}, shift::Number)
    for (i,bus) in network["bus"]
        bus["va"] = bus["va"] + shift
    end
end


"add start values to a data model"
function c1_set_start_values!(network::Dict{String,<:Any}; branch_flow=false)
    for (i,bus) in network["bus"]
        bus["va_start"] = bus["va"]
        bus["vm_start"] = bus["vm"]

        bus["vr_start"] = bus["vm"]*cos(bus["va"])
        bus["vi_start"] = bus["vm"]*sin(bus["va"])

        if haskey(bus, "vm_offset") && !isnan(bus["vm_offset"])
            bus["vm_offset_start"] = bus["vm_offset"]
        end
    end

    for (i,gen) in network["gen"]
        gen["pg_start"] = gen["pg"]
        gen["qg_start"] = gen["qg"]
    end

    for (i,shunt) in network["shunt"]
        shunt["bs_start"] = shunt["bs"]
    end

    if branch_flow
        for (i,branch) in network["branch"]
            branch["pf_start"] = branch["pf"]
            branch["qf_start"] = branch["qf"]

            branch["pt_start"] = branch["pt"]
            branch["qt_start"] = branch["qt"]
        end
    end
end


function update_active_power_data!(network::Dict{String,<:Any}, data::Dict{String,<:Any}; branch_flow=false)
    for (i,bus) in data["bus"]
        nw_bus = network["bus"][i]
        nw_bus["va"] = bus["va"]
    end

    for (i,gen) in data["gen"]
        nw_gen = network["gen"][i]
        nw_gen["pg"] = gen["pg"]
    end

    if branch_flow
        for (i,branch) in data["branch"]
            nw_branch = network["branch"][i]
            nw_branch["pf"] = branch["pf"]
            nw_branch["pt"] = branch["pt"]
        end
    end
end


function c1_extract_solution(network::Dict{String,<:Any}; branch_flow=false)
    sol = Dict{String,Any}()

    sol["bus"] = Dict{String,Any}()
    for (i,bus) in network["bus"]
        bus_dict = Dict{String,Any}()
        bus_dict["va"] = get(bus, "va", 0.0)
        bus_dict["vm"] = get(bus, "vm", 1.0)
        sol["bus"][i] = bus_dict
    end

    sol["shunt"] = Dict{String,Any}()
    for (i,shunt) in network["shunt"]
        shunt_dict = Dict{String,Any}()
        shunt_dict["gs"] = get(shunt, "gs", 0.0)
        shunt_dict["bs"] = get(shunt, "bs", 0.0)
        sol["shunt"][i] = shunt_dict
    end

    sol["gen"] = Dict{String,Any}()
    for (i,gen) in network["gen"]
        gen_dict = Dict{String,Any}()
        gen_dict["pg"] = get(gen, "pg", 0.0)
        gen_dict["qg"] = get(gen, "qg", 0.0)
        sol["gen"][i] = gen_dict
    end

    if branch_flow
        sol["branch"] = Dict{String,Any}()
        for (i,branch) in network["branch"]
            branch_dict = Dict{String,Any}()
            branch_dict["pf"] = get(branch, "pf", 0.0)
            branch_dict["qf"] = get(branch, "qf", 0.0)
            branch_dict["pt"] = get(branch, "pt", 0.0)
            branch_dict["qt"] = get(branch, "qt", 0.0)
            sol["branch"][i] = branch_dict
        end
    end

    return sol
end



##### Solution Analysis #####


"assumes a vaild ac solution is included in the data and computes the branch flow values"
function calc_c1_branch_flow_ac(data::Dict{String,<:Any})
    @assert("per_unit" in keys(data) && data["per_unit"])
    @assert(!haskey(data, "conductors"))

    if _IM.ismultinetwork(data)
        nws = Dict{String,Any}()
        for (i,nw_data) in data["nw"]
            nws[i] = _calc_c1_branch_flow_ac(nw_data)
        end
        return Dict{String,Any}(
            "nw" => nws,
            "per_unit" => data["per_unit"],
            "baseMVA" => data["baseMVA"]
        )
    else
        flows = _calc_c1_branch_flow_ac(data)
        flows["per_unit"] = data["per_unit"]
        flows["baseMVA"] = data["baseMVA"]
        return flows
    end
end


"helper function for calc_branch_flow_ac"
function _calc_c1_branch_flow_ac(data::Dict{String,<:Any})
    vm = Dict(bus["index"] => bus["vm"] for (i,bus) in data["bus"])
    va = Dict(bus["index"] => bus["va"] for (i,bus) in data["bus"])

    flows = Dict{String,Any}()
    for (i,branch) in data["branch"]
        if branch["br_status"] != 0
            f_bus = branch["f_bus"]
            t_bus = branch["t_bus"]

            g, b = _PM.calc_branch_y(branch)
            tr, ti = _PM.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]
            g_to = branch["g_to"]
            b_to = branch["b_to"]

            tm = branch["tap"]

            vm_fr = vm[f_bus]
            vm_to = vm[t_bus]
            va_fr = va[f_bus]
            va_to = va[t_bus]

            if branch["transformer"]
                p_fr =  (g/tm^2+g_fr)*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))
                q_fr = -(b/tm^2+b_fr)*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))
            else
                p_fr =  (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))
                q_fr = -(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))
            end

            p_to =  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))
            q_to = -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))
        else
            p_fr = NaN
            q_fr = NaN

            p_to = NaN
            q_to = NaN
        end

        flows[i] = Dict(
            "pf" => p_fr,
            "qf" => q_fr,
            "pt" => p_to,
            "qt" => q_to
        )
    end

    return Dict{String,Any}("branch" => flows)
end

function calc_c1_power_balance_deltas!(network::Dict{String,<:Any})
    flows = calc_c1_branch_flow_ac(network)
    _PM.update_data!(network, flows)

    balance = _PM.calc_power_balance(network)
    _PM.update_data!(network, balance)

    p_delta_abs = [abs(bus["p_delta"]) for (i,bus) in network["bus"] if bus["bus_type"] != 4]
    q_delta_abs = [abs(bus["q_delta"]) for (i,bus) in network["bus"] if bus["bus_type"] != 4]

    return (
        p_delta_abs_max = maximum(p_delta_abs),
        p_delta_abs_mean = mean(p_delta_abs),
        q_delta_abs_max = maximum(q_delta_abs),
        q_delta_abs_mean = mean(q_delta_abs),
    )
end



function calc_c1_violations(network::Dict{String,<:Any}, solution::Dict{String,<:Any}; vm_digits=3, rate_key="rate_c")
    vm_vio = 0.0
    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            bus_sol = solution["bus"][i]

            # helps to account for minor errors in equality constraints
            sol_val = round(bus_sol["vm"], digits=vm_digits)

            #vio_flag = false
            if sol_val < bus["vmin"]
                vm_vio += bus["vmin"] - sol_val
                #vio_flag = true
            end
            if sol_val > bus["vmax"]
                vm_vio += sol_val - bus["vmax"]
                #vio_flag = true
            end
            #if vio_flag
            #    info(_LOGGER, "$(i): $(bus["vmin"]) - $(sol_val) - $(bus["vmax"])")
            #end
        end
    end

    pg_vio = 0.0
    qg_vio = 0.0
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            gen_sol = solution["gen"][i]

            if gen_sol["pg"] < gen["pmin"]
                pg_vio += gen["pmin"] - gen_sol["pg"]
            end
            if gen_sol["pg"] > gen["pmax"]
                pg_vio += gen_sol["pg"] - gen["pmax"]
            end

            if gen_sol["qg"] < gen["qmin"]
                qg_vio += gen["qmin"] - gen_sol["qg"]
            end
            if gen_sol["qg"] > gen["qmax"]
                qg_vio += gen_sol["qg"] - gen["qmax"]
            end
        end
    end


    sm_vio = NaN
    if haskey(solution, "branch")
        sm_vio = 0.0
        for (i,branch) in network["branch"]
            if branch["br_status"] != 0
                branch_sol = solution["branch"][i]

                s_fr = abs(branch_sol["pf"])
                s_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                # note true model is rate_c
                #vio_flag = false
                rating = branch[rate_key]

                if s_fr > rating
                    sm_vio += s_fr - rating
                    #vio_flag = true
                end
                if s_to > rating
                    sm_vio += s_to - rating
                    #vio_flag = true
                end
                #if vio_flag
                #    info(_LOGGER, "$(i), $(branch["f_bus"]), $(branch["t_bus"]): $(s_fr) / $(s_to) <= $(branch["rate_c"])")
                #end
            end
        end
    end

    return (vm=vm_vio, pg=pg_vio, qg=qg_vio, sm=sm_vio)
end


"returns a sorted list of branch flow violations"
function branch_c1_violations_sorted(network::Dict{String,<:Any}, solution::Dict{String,<:Any}; rate_key="rate_c")
    branch_violations = []

    if haskey(solution, "branch")
        for (i,branch) in network["branch"]
            if branch["br_status"] != 0
                branch_sol = solution["branch"][i]

                s_fr = abs(branch_sol["pf"])
                s_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                sm_vio = 0.0

                rating = branch[rate_key]
                if s_fr > rating
                    sm_vio = s_fr - rating
                end
                if s_to > rating && s_to - rating > sm_vio
                    sm_vio = s_to - rating
                end

                if sm_vio > 0.0
                    push!(branch_violations, (branch_id=branch["index"], sm_vio=sm_vio))
                end
            end
        end
    end

    sort!(branch_violations, by=(x) -> -x.sm_vio)

    return branch_violations
end







"updates tap changing transformers to maintain nominal voltage values"
function normalize_tm_values!(data)
    bus_voltage_range = Dict{Int,Any}()
    for (i,bus) in data["bus"]
        bus_voltage_range[bus["index"]] = (vmin=bus["vmin"], vmax=bus["vmax"])
    end

    for (i,branch) in data["branch"]
        if branch["transformer"] && abs(branch["control_mode"]) == 1
            vm_target_fr = sum(bus_voltage_range[branch["f_bus"]])/2.0
            vm_target_to = sum(bus_voltage_range[branch["t_bus"]])/2.0

            tap_target = vm_target_fr/vm_target_to

            if tap_target < branch["tm_min"]
                info(_LOGGER, "transformer branch $(i) unable to hit tap target $(tap_target) -> $(branch["tm_min"])")
                tap_target = branch["tm_min"]
            end

            if tap_target > branch["tm_max"]
                info(_LOGGER, "transformer branch $(i) unable to hit tap target $(tap_target) -> $(branch["tm_max"])")
                tap_target = branch["tm_max"]
            end

            tm_step, tm_value = closest_discrete_value(tap_target, branch["tm_min"], branch["tm_max"], branch["tm_steps"])

            if tm_step != branch["tm_step"]
                info(_LOGGER, "transformer branch $(i) tap, $(branch["tap"])/$(branch["tm_step"])  -> $(tm_value)/$(tm_step)")
                branch["tap"] = tm_value
                branch["tm_step"] = tm_step

                if branch["tab1"] != 0
                    ict = data["impedance_correction_table"][branch["tab1"]]
                    correct_transformer_impedance!(branch, ict, branch["tap"])
                end
            end
        end
    end
end

"converts a continuous value to its nearest integer value"
function closest_discrete_value(value::Real, range_lb::Real, range_ub::Real, steps::Int)
    @assert(range_lb <= value && value <= range_ub)
    @assert(steps >= 0.0)

    mid_point = (range_lb + range_ub)/2.0
    target_value = value - mid_point
    step_size = (range_ub - range_lb)/(steps-1)
    target_step = round(Int, target_value/step_size)

    closest_value = mid_point + target_step*step_size

    return target_step, closest_value
end

""
function set_va_start_values!(data)
    for (i,bus) in data["bus"]
        bus["va_start"] = bus["va"]
    end
end


"update branch impedance values based on the given impedance correction points and tap setting"
function correct_transformer_impedance!(branch, impedance_correction_points, tap_setting::Real)
    settings = comp_transformer_impedance!(branch, impedance_correction_points, tap_setting)

    branch["br_r"] = settings["br_r"]
    branch["br_x"] = settings["br_x"]
end

"compute the branch impedance values based on the given impedance correction points and tap setting"
function comp_transformer_impedance!(branch, impedance_correction_points, tap_setting::Real)

    if tap_setting < impedance_correction_points[1].t || tap_setting > impedance_correction_points[end].t
        error(_LOGGER, "tap setting out of bounds of the impedance correction data")
    end

    icp1 = impedance_correction_points[1]
    icp2 = impedance_correction_points[2]

    for i in 1:(length(impedance_correction_points)-1)
        icp1 = impedance_correction_points[i]
        icp2 = impedance_correction_points[i+1]
        if icp1.t <= tap_setting && tap_setting <= icp2.t
            #println("$(icp1.t) - $(tap_setting) - $(icp2.t)  OK!")
            break
        else
            #println("$(icp1.t) - $(tap_setting) - $(icp2.t)")
        end
    end

    x1 = icp1.t
    y1 = icp1.f
    x2 = icp2.t
    y2 = icp2.f

    m = (y2 - y1)/(x2 - x1)
    b = y1 - m * x1

    scaling = m*tap_setting + b

    #println(icp1, " ", icp2, " ", tap_setting, " ", scaling)
    return Dict(
        "br_r" => branch["br_r_nominal"]*scaling,
        "br_x" => branch["br_x_nominal"]*scaling
    )
end
