##### Shared PowerModels Data Transformation #####


"""
transforms files from ARPA-e GOC Challenge 1 data format in to the PowerModels
data format.  This consists of taking the data from multiple data structures
and putting into a network data dictionary.
"""
function build_pm_model(goc_data)
    scenario = goc_data.scenario
    network = goc_data.network

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
    for dispatch_tbl in goc_data.cost["disptbl"]
        dispatch_tbl_lookup[dispatch_tbl["ctbl"]] = dispatch_tbl
    end

    cost_tbl_lookup = Dict()
    for cost_tbl in goc_data.cost["ctbl"]
        cost_tbl_lookup[cost_tbl["ltbl"]] = cost_tbl
    end

    gen_cost_models = Dict()
    for gen_dispatch in goc_data.cost["gen"]
        gen_id = (gen_dispatch["bus"], strip(replace(gen_dispatch["genid"], "'" => "")))
        dispatch_tbl = dispatch_tbl_lookup[gen_dispatch["disptbl"]]
        cost_tbl = cost_tbl_lookup[dispatch_tbl["ctbl"]]

        gen_cost_models[gen_id] = cost_tbl
    end

    if length(gen_cost_models) != length(network["gen"])
        error(LOGGER, "cost model data missing, network has $(length(network["gen"])) generators, the cost model has $(length(gen_cost_models)) generators")
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

    if length(goc_data.response) != length(network["gen"])
        error(LOGGER, "generator response model data missing, network has $(length(network["gen"])) generators, the response model has $(length(goc_data.response)) generators")
    end

    for gen_response in goc_data.response
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
            error(LOGGER, "unrecognized contingency component type $(cont["component"]) at contingency $(i)")
        end
    end

    network["branch_contingencies"] = branch_ids
    network["gen_contingencies"] = generator_ids

    network["branch_contingencies_active"] = []
    network["gen_contingencies_active"] = []



    ##### Fix Broken Data #####

    PowerModels.correct_cost_functions!(network)

    # FYI, this breaks output API
    #PowerModels.propagate_topology_status!(network)

    for (i,shunt) in network["shunt"]
        # test checks if a "switched shunt" in the orginal data model
        if shunt["dispatchable"]
            if shunt["bs"] < shunt["bmin"]
                warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmin"])")
                shunt["bs"] = shunt["bmin"]
            end
            if shunt["bs"] > shunt["bmax"]
                warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmax"])")
                shunt["bs"] = shunt["bmax"]
            end
        end
    end

    return network
end


function build_pm_opf_model(goc_data)
    scenario = goc_data.scenario
    network = goc_data.network

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
    for dispatch_tbl in goc_data.cost["disptbl"]
        dispatch_tbl_lookup[dispatch_tbl["ctbl"]] = dispatch_tbl
    end

    cost_tbl_lookup = Dict()
    for cost_tbl in goc_data.cost["ctbl"]
        cost_tbl_lookup[cost_tbl["ltbl"]] = cost_tbl
    end

    gen_cost_models = Dict()
    for gen_dispatch in goc_data.cost["gen"]
        gen_id = (gen_dispatch["bus"], strip(gen_dispatch["genid"]))
        dispatch_tbl = dispatch_tbl_lookup[gen_dispatch["disptbl"]]
        cost_tbl = cost_tbl_lookup[dispatch_tbl["ctbl"]]

        gen_cost_models[gen_id] = cost_tbl
    end

    if length(gen_cost_models) != length(network["gen"])
        error(LOGGER, "cost model data missing, network has $(length(network["gen"])) generators, the cost model has $(length(gen_cost_models)) generators")
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

    PowerModels.correct_cost_functions!(network)

    # FYI, this breaks output API
    #PowerModels.propagate_topology_status!(network)

    for (i,shunt) in network["shunt"]
        # test checks if a "switched shunt" in the orginal data model
        if shunt["dispatchable"]
            if shunt["bs"] < shunt["bmin"]
                warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmin"])")
                shunt["bs"] = shunt["bmin"]
            end
            if shunt["bs"] > shunt["bmax"]
                warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmax"])")
                shunt["bs"] = shunt["bmax"]
            end
        end
    end


    return network
end



function build_scopf_multinetwork(network)
    if InfrastructureModels.ismultinetwork(network)
        error(LOGGER, "build scopf can only be used on single networks")
    end

    contingencies = length(network["gen_contingencies"]) + length(network["branch_contingencies"])

    info(LOGGER, "building scopf multi-network with $(contingencies) networks")

    mn_data = PowerModels.replicate(network, contingencies)
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

        #=
        # TODO restore this
        for (i,branch) in network["branch"]
            if haskey(branch, "rate_c")
                branch["rate_a"] = branch["rate_c"]
            end
        end
        =#
    end

    network_id = 1
    for cont in base_network["gen_contingencies"]
        cont_nw = mn_data["nw"]["$(network_id)"]
        cont_nw["name"] = cont.label
        cont_nw["gen"]["$(cont.idx)"]["gen_status"] = 0
        network_id += 1
    end
    for cont in base_network["branch_contingencies"]
        cont_nw = mn_data["nw"]["$(network_id)"]
        cont_nw["name"] = cont.label
        cont_nw["branch"]["$(cont.idx)"]["br_status"] = 0
        network_id += 1
    end

    return mn_data
end



function read_solution1(network; output_dir="", state_file="solution1.txt")
    if length(output_dir) > 0
        solution1_path = joinpath(output_dir, state_file)
    else
        solution1_path = state_file
    end

    return build_pm_solution(network, solution1_path)
end

function build_pm_solution(network, goc_sol_file::String)
    info(LOGGER, "loading solution file: $(goc_sol_file)")
    goc_sol = parse_solution1_file(goc_sol_file)

    info(LOGGER, "converting GOC solution to PowerModels solution")
    pm_sol = build_pm_solution(network, goc_sol)

    return pm_sol
end

function build_pm_solution(network, goc_sol)
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



##### Shared Data Transformations #####

function tighten_constraints!(pm_network)
    for (i,bus) in pm_network["bus"]
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

    for (i,gen) in pm_network["gen"]
        if gen["pmax"] > 0 && gen["pmax"]*0.9 > gen["pmin"]
            gen["pmax_target"] = gen["pmax"]*0.9
        else
            gen["pmax_target"] = gen["pmax"]
        end
    end

    for (i,branch) in pm_network["branch"]
        # TODO restore this
        #if isapprox(branch["rate_a"], branch["rate_c"])
            branch["rate_a"] = branch["rate_a"]*0.80
        #end
    end
end



function deactivate_rate_a!(pm_network)
    pm_network["active_rates"] = Int[]
    for (i,branch) in pm_network["branch"]
        branch["rate_a_inactive"] = branch["rate_a"]
        delete!(branch, "rate_a")
    end
end

function activate_rate_a!(pm_network)
    if haskey(pm_network, "active_rates")
        delete!(pm_network, "active_rates")
    end

    for (i,branch) in pm_network["branch"]
        if haskey(branch, "rate_a_inactive")
            branch["rate_a"] = branch["rate_a_inactive"]
            delete!(branch, "rate_a_inactive")
        end
    end
end

function activate_rate_a_violations!(pm_network)
    ac_flows = PowerModels.calc_branch_flow_ac(pm_network)
    for (i,branch) in pm_network["branch"]
        branch["pf_start"] = ac_flows["branch"][i]["pf"]
        branch["qf_start"] = ac_flows["branch"][i]["qf"]

        branch["pt_start"] = ac_flows["branch"][i]["pt"]
        branch["qt_start"] = ac_flows["branch"][i]["qt"]
    end

    line_flow_vio = false
    for (i,branch) in pm_network["branch"]
        if !haskey(branch, "rate_a")
            if (ac_flows["branch"][i]["pf"]^2 + ac_flows["branch"][i]["qf"]^2 > branch["rate_a_inactive"]^2 ||
                ac_flows["branch"][i]["pt"]^2 + ac_flows["branch"][i]["qt"]^2 > branch["rate_a_inactive"]^2)
                info(LOGGER, "add rate_a flow limit on branch $(i) $(branch["source_id"])")
                #branch["rate_a"] = branch["rate_a_inactive"] - max(abs(ac_flows["branch"][i]["qf"]), abs(ac_flows["branch"][i]["qt"]))
                branch["rate_a"] = branch["rate_a_inactive"]
                push!(pm_network["active_rates"], branch["index"])
                line_flow_vio = true
            end
        else
            sm_fr = sqrt(ac_flows["branch"][i]["pf"]^2 + ac_flows["branch"][i]["qf"]^2)
            sm_to = sqrt(ac_flows["branch"][i]["pf"]^2 + ac_flows["branch"][i]["qf"]^2)
            vio = max(0.0, sm_fr - branch["rate_a"], sm_to - branch["rate_a"])
            if vio > 0.01
                warn(LOGGER, "add rate_a flow limit violations $(vio) on branch $(i) $(branch["source_id"])")
            end
        end
    end

    return line_flow_vio
end



function gens_by_bus(network)
    bus_gens = Dict(i => Any[] for (i,bus) in network["bus"])
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            push!(bus_gens["$(gen["gen_bus"])"], gen)
        end
    end
    return bus_gens
end


"build a static ordering of all contigencies"
function contingency_order(pm_network)
    gen_cont_order = sort(pm_network["gen_contingencies"], by=(x) -> x.label)
    branch_cont_order = sort(pm_network["branch_contingencies"], by=(x) -> x.label)

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



"""
assumes there is one reference bus and one connected component and adjusts voltage
angles to be centered around zero at the reference bus.
"""
function correct_voltage_angles!(network)
    ref_bus = -1
    for (i,bus) in network["bus"]
        if bus["bus_type"] == 3
            @assert ref_bus == -1
            ref_bus = bus
        end
    end

    if !isapprox(ref_bus["va"], 0.0, atol=1e-8)
        warn(LOGGER, "shifting voltage angles by $(-ref_bus["va"]) to set reference bus to 0.0")
        shift_voltage_anlges!(network, -ref_bus["va"])
    end
end


"shift networks voltage angles by a specified amount"
function shift_voltage_anlges!(network, shift::Number)
    for (i,bus) in network["bus"]
        bus["va"] = bus["va"] + shift
    end
end


##### Shared PowerModels Extensions #####



"generates variables for both `active` and `reactive` bus deltas"
function variable_bus_delta_abs(pm::AbstractPowerModel; kwargs...)
    variable_active_delta_abs(pm; kwargs...)
    variable_reactive_delta_abs(pm; kwargs...)
end


""
function variable_active_delta_abs(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    if bounded
         var(pm, nw, cnd)[:p_delta_abs] = @variable(pm.model,
            [i in ids(pm, :bus)], base_name="$(nw)_$(cnd)_p_delta_abs",
            upper_bound = 0.5,
            lower_bound = 0,
            start = 0.0
        )
    else
         var(pm, nw, cnd)[:p_delta_abs] = @variable(pm.model,
            [i in ids(pm, :bus)], base_name="$(nw)_$(cnd)_p_delta_abs",
            start = 0.0
        )
    end
end

""
function variable_reactive_delta_abs(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    if bounded
         var(pm, nw, cnd)[:q_delta_abs] = @variable(pm.model,
            [i in ids(pm, :bus)], base_name="$(nw)_$(cnd)_q_delta_abs",
            upper_bound = 0.5,
            lower_bound = 0,
            start = 0.0
        )
    else
         var(pm, nw, cnd)[:q_delta_abs] = @variable(pm.model,
            [i in ids(pm, :bus)], base_name="$(nw)_$(cnd)_q_delta_abs",
            start = 0.0
        )
    end
end


""
function variable_reactive_shunt(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    var(pm, nw, cnd)[:bsh] = @variable(pm.model,
        [i in ids(pm, nw, :shunt_var)], base_name="$(nw)_$(cnd)_bsh",
        upper_bound = ref(pm, nw, :shunt, i, "bmax", cnd),
        lower_bound = ref(pm, nw, :shunt, i, "bmin", cnd),
        start = comp_start_value(ref(pm, nw, :shunt, i), "bsh_start", cnd)
    )
end

""
function variable_reactive_shunt(pm::AbstractWModels; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    var(pm, nw, cnd)[:bsh] = @variable(pm.model,
        [i in ids(pm, nw, :shunt_var)], base_name="$(nw)_$(cnd)_bsh",
        upper_bound = ref(pm, nw, :shunt, i, "bmax", cnd),
        lower_bound = ref(pm, nw, :shunt, i, "bmin", cnd),
        start = comp_start_value(ref(pm, nw, :shunt, i), "bsh_start", cnd)
    )

    var(pm, nw, cnd)[:wbsh] = @variable(pm.model,
        [i in ids(pm, nw, :shunt_var)], base_name="$(nw)_$(cnd)_wbsh",
        start = 0.0
    )
end


""
function constraint_power_balance_shunt_dispatch(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(con(pm, nw, cnd), :kcl_p)
        con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(con(pm, nw, cnd), :kcl_q)
        con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts_const = ref(pm, :bus_shunts_const, i)
    bus_shunts_var = ref(pm, :bus_shunts_var, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs_const = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts_const)
    bus_bs_const = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts_const)

    constraint_power_balance_shunt_dispatch(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
end

""
function constraint_power_balance_shunt_dispatch(pm::AbstractACPModel, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vm = var(pm, n, c, :vm, i)
    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    bsh = var(pm, n, c, :bsh)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    con(pm, n, c, :kcl_p)[i] = JuMP.@NLconstraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm^2)
    con(pm, n, c, :kcl_q)[i] = JuMP.@NLconstraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm^2 + sum(bsh[s]*vm^2 for s in bus_shunts_var))
end

""
function constraint_power_balance_shunt_dispatch(pm::AbstractACRModel, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vi = var(pm, n, c, :vi, i)
    vr = var(pm, n, c, :vr, i)
    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    bsh = var(pm, n, c, :bsh)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    # possibly can save 2x in function eval, but no the dominant runtime at this moment
    #vm_sqr = @variable(pm.model, start=1.0, base_name="$(0)_vm_sqr_$(i)")

    #JuMP.@constraint(pm.model, vm_sqr == vi^2 + vr^2)
    #con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm_sqr)
    #con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm_sqr + sum(bsh[s]*vm_sqr for s in bus_shunts_var))

    con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*(vi^2 + vr^2))
    con(pm, n, c, :kcl_q)[i] = JuMP.@NLconstraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*(vi^2 + vr^2) + sum(bsh[s]*(vi^2 + vr^2) for s in bus_shunts_var))
end

""
function constraint_power_balance_shunt_dispatch(pm::AbstractWRModels, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    w = var(pm, n, c, :w, i)
    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    bsh = var(pm, n, c, :bsh)
    wbsh = var(pm, n, c, :wbsh)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w)
    con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbsh[s] for s in bus_shunts_var))

    for s in bus_shunts_var
        InfrastructureModels.relaxation_product(pm.model, w, bsh[s], wbsh[s])
    end
end

""
function constraint_power_balance_shunt_dispatch(pm::AbstractActivePowerModel, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    p = var(pm, n, c, :p)
    pg = var(pm, n, c, :pg)
    p_dc = var(pm, n, c, :p_dc)

    con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*1.0)
end


""
function constraint_power_balance_shunt_dispatch_soft(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(con(pm, nw, cnd), :kcl_p)
        con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(con(pm, nw, cnd), :kcl_q)
        con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts_const = ref(pm, :bus_shunts_const, i)
    bus_shunts_var = ref(pm, :bus_shunts_var, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs_const = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts_const)
    bus_bs_const = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts_const)

    constraint_power_balance_shunt_dispatch_soft(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
end

""
function constraint_power_balance_shunt_dispatch_soft(pm::AbstractACPModel, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vm = var(pm, n, c, :vm, i)
    p_delta_abs = var(pm, n, c, :p_delta_abs, i)
    q_delta_abs = var(pm, n, c, :q_delta_abs, i)

    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    bsh = var(pm, n, c, :bsh)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    p_delta = @NLexpression(pm.model, - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm^2)
    q_delta = @NLexpression(pm.model, - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm^2 + sum(bsh[s]*vm^2 for s in bus_shunts_var))

    @NLconstraint(pm.model,  p_delta_abs >= p_delta)
    @NLconstraint(pm.model, -p_delta_abs <= p_delta)

    @NLconstraint(pm.model,  q_delta_abs >= q_delta)
    @NLconstraint(pm.model, -q_delta_abs <= q_delta)
end

""
function constraint_power_balance_shunt_dispatch_soft(pm::AbstractWRModels, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    w = var(pm, n, c, :w, i)
    p_delta_abs = var(pm, n, c, :p_delta_abs, i)
    q_delta_abs = var(pm, n, c, :q_delta_abs, i)

    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    bsh = var(pm, n, c, :bsh)
    wbsh = var(pm, n, c, :wbsh)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    #p_delta = - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w
    #q_delta = - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbsh[s] for s in bus_shunts_var)

    @constraint(pm.model,  p_delta_abs >= - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w)
    @constraint(pm.model, -p_delta_abs <= - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*w)

    @constraint(pm.model,  q_delta_abs >= - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbsh[s] for s in bus_shunts_var))
    @constraint(pm.model, -q_delta_abs <= - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*w + sum(wbsh[s] for s in bus_shunts_var))

    for s in bus_shunts_var
        InfrastructureModels.relaxation_product(pm.model, w, bsh[s], wbsh[s])
    end
end



""
function constraint_ohms_yt_from_goc(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"][cnd]
    b_fr = branch["b_fr"][cnd]
    tm = branch["tap"][cnd]

    if branch["transformer"]
        constraint_ohms_yt_from_goc(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g[cnd,cnd], b[cnd,cnd], g_fr, b_fr, tr[cnd], ti[cnd], tm)
    else
        PowerModels.constraint_ohms_yt_from(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g[cnd,cnd], b[cnd,cnd], g_fr, b_fr, tr[cnd], ti[cnd], tm)
    end
end

function constraint_ohms_yt_from_goc(pm::AbstractACPModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n, c,  :p, f_idx)
    q_fr  = var(pm, n, c,  :q, f_idx)
    vm_fr = var(pm, n, c, :vm, f_bus)
    vm_to = var(pm, n, c, :vm, t_bus)
    va_fr = var(pm, n, c, :va, f_bus)
    va_to = var(pm, n, c, :va, t_bus)

    JuMP.@NLconstraint(pm.model, p_fr ==  (g/tm^2+g_fr)*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
    JuMP.@NLconstraint(pm.model, q_fr == -(b/tm^2+b_fr)*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to)) )
end

function constraint_ohms_yt_from_goc(pm::AbstractACRModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, c, :p, f_idx)
    q_fr = var(pm, n, c, :q, f_idx)
    vr_fr = var(pm, n, c, :vr, f_bus)
    vr_to = var(pm, n, c, :vr, t_bus)
    vi_fr = var(pm, n, c, :vi, f_bus)
    vi_to = var(pm, n, c, :vi, t_bus)

    JuMP.@constraint(pm.model, p_fr ==  (g/tm^2+g_fr)*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
    JuMP.@constraint(pm.model, q_fr == -(b/tm^2+b_fr)*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
end

function constraint_ohms_yt_from_goc(pm::AbstractWRModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr = var(pm, n, c, :p, f_idx)
    q_fr = var(pm, n, c, :q, f_idx)
    w_fr = var(pm, n, c, :w, f_bus)
    wr   = var(pm, n, c, :wr, (f_bus, t_bus))
    wi   = var(pm, n, c, :wi, (f_bus, t_bus))

    JuMP.@constraint(pm.model, p_fr ==  (g/tm^2+g_fr)*w_fr + (-g*tr+b*ti)/tm^2*wr + (-b*tr-g*ti)/tm^2*wi )
    JuMP.@constraint(pm.model, q_fr == -(b/tm^2+b_fr)*w_fr - (-b*tr-g*ti)/tm^2*wr + (-g*tr+b*ti)/tm^2*wi )
end

function constraint_ohms_yt_from_goc(pm::AbstractDCPModel, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    p_fr  = var(pm, n, c,  :p, f_idx)
    va_fr = var(pm, n, c, :va, f_bus)
    va_to = var(pm, n, c, :va, t_bus)

    JuMP.@constraint(pm.model, p_fr == -b*(va_fr - va_to))
    # omit reactive constraint
end


"adds pg_cost variables and constraints"
function objective_variable_pg_cost(pm::AbstractPowerModel)
    for (n, nw_ref) in nws(pm)
        gen_lines = calc_cost_pwl_lines(nw_ref[:gen])
        pg_cost_start = Dict{Int64,Float64}()
        for (i, gen) in nw_ref[:gen]
            pg_value = sum(JuMP.start_value(var(pm, n, c, :pg, i)) for c in conductor_ids(pm, n))
            pg_cost_value = -Inf
            for line in gen_lines[i]
                pg_cost_value = max(pg_cost_value, line.slope*pg_value + line.intercept)
            end
            pg_cost_start[i] = pg_cost_value
        end

        #println(pg_cost_start)
        #println(sum(values(pg_cost_start)))

        pg_cost = var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in ids(pm, n, :gen)], base_name="$(n)_pg_cost",
            start=pg_cost_start[i]
        )

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            for line in gen_lines[i]
                JuMP.@constraint(pm.model, pg_cost[i] >= line.slope*sum(var(pm, n, c, :pg, i) for c in conductor_ids(pm, n)) + line.intercept)
            end
        end
    end
end


"add start values to a data model"
function set_start_values!(network::Dict{String,<:Any}; branch_flow=false)
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
        shunt["bsh_start"] = shunt["bs"]
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


function ref_add_goc!(pm::AbstractPowerModel)
    if InfrastructureModels.ismultinetwork(pm.data)
        nws_data = pm.data["nw"]
    else
        nws_data = Dict("0" => pm.data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        ref = pm.ref[:nw][nw_id]

        ref[:shunt_const] = Dict(x for x in ref[:shunt] if (!haskey(x.second, "dispatchable") || !x.second["dispatchable"]))
        ref[:shunt_var] = Dict(x for x in ref[:shunt] if (haskey(x.second, "dispatchable") && x.second["dispatchable"]))

        bus_shunts_const = Dict((i, []) for (i,bus) in ref[:bus])
        for (i,shunt) in ref[:shunt_const]
            push!(bus_shunts_const[shunt["shunt_bus"]], i)
        end
        ref[:bus_shunts_const] = bus_shunts_const

        bus_shunts_var = Dict((i, []) for (i,bus) in ref[:bus])
        for (i,shunt) in ref[:shunt_var]
            push!(bus_shunts_var[shunt["shunt_bus"]], i)
        end
        ref[:bus_shunts_var] = bus_shunts_var
    end
end


function add_setpoint_dispatchable(
    sol,
    pm::AbstractPowerModel,
    dict_name,
    param_name,
    variable_symbol;
    index_name = "index",
    default_value = (item) -> NaN,
    scale = (x,item,cnd) -> x,
    extract_var = (var,idx,item) -> var[idx],
    sol_dict = get(sol, dict_name, Dict{String,Any}()),
    conductorless = false,
    dispatchable_check = false
)

    if InfrastructureModels.ismultinetwork(pm.data)
        data_dict = pm.data["nw"]["$(pm.cnw)"][dict_name]
    else
        data_dict = pm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end
    for (i,item) in data_dict
        if dispatchable_check && (!haskey(item, "dispatchable") || !item["dispatchable"])
            continue
        end

        idx = Int(item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())

        if conductorless
            sol_item[param_name] = default_value(item)
            try
                variable = extract_var(var(pm, pm.cnw, variable_symbol), idx, item)
                sol_item[param_name] = scale(JuMP.value(variable), item, 1)
            catch
            end
        else
            num_conductors = length(conductor_ids(pm))
            cnd_idx = 1
            sol_item[param_name] = MultiConductorVector{Real}([default_value(item) for i in 1:num_conductors])
            for conductor in conductor_ids(pm)
                try
                    variable = extract_var(var(pm, variable_symbol, cnd=conductor), idx, item)
                    sol_item[param_name][cnd_idx] = scale(JuMP.value(variable), item, conductor)
                catch
                end
                cnd_idx += 1
            end
        end

        # remove MultiConductorValue, if it was not a ismulticonductor network
        if !ismulticonductor(pm)
            sol_item[param_name] = sol_item[param_name][1]
        end
    end
end


function update_active_power_data!(network, data; branch_flow=false)
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


function extract_solution(network; branch_flow=false)
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



##### GOC Solution Analysis #####

function compute_violations(network, solution; vm_digits=3)

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
            #    info(LOGGER, "$(i): $(bus["vmin"]) - $(sol_val) - $(bus["vmax"])")
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
                s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)

                # note true model is rate_c
                #vio_flag = false
                if s_fr > branch["rate_a"]
                    sm_vio += s_fr - branch["rate_a"]
                    #vio_flag = true
                end
                if s_to > branch["rate_a"]
                    sm_vio += s_to - branch["rate_a"]
                    #vio_flag = true
                end

                #=
                if s_fr > branch["rate_c"]
                    sm_vio += s_fr - branch["rate_c"]
                    #vio_flag = true
                end
                if s_to > branch["rate_c"]
                    sm_vio += s_to - branch["rate_c"]
                    #vio_flag = true
                end
                =#
                #if vio_flag
                #    info(LOGGER, "$(i), $(branch["f_bus"]), $(branch["t_bus"]): $(s_fr) / $(s_to) <= $(branch["rate_c"])")
                #end
            end
        end
    end

    return (vm=vm_vio, pg=pg_vio, qg=qg_vio, sm=sm_vio)
end



function compute_violations_ratec(network, solution; vm_digits=3)

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
            #    info(LOGGER, "$(i): $(bus["vmin"]) - $(sol_val) - $(bus["vmax"])")
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
                s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)

                # note true model is rate_c
                #vio_flag = false
                if s_fr > branch["rate_c"]
                    sm_vio += s_fr - branch["rate_c"]
                    #vio_flag = true
                end
                if s_to > branch["rate_c"]
                    sm_vio += s_to - branch["rate_c"]
                    #vio_flag = true
                end
                #if vio_flag
                #    info(LOGGER, "$(i), $(branch["f_bus"]), $(branch["t_bus"]): $(s_fr) / $(s_to) <= $(branch["rate_c"])")
                #end
            end
        end
    end

    return (vm=vm_vio, pg=pg_vio, qg=qg_vio, sm=sm_vio)
end

"returns a sorted list of branch flow violations"
function branch_violations_sorted(network, solution)
    branch_violations = []

    if haskey(solution, "branch")
        for (i,branch) in network["branch"]
            if branch["br_status"] != 0
                branch_sol = solution["branch"][i]
                s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)

                sm_vio = 0.0
                # TODO update to rate_c
                if s_fr > branch["rate_a"]
                    sm_vio = s_fr - branch["rate_a"]
                end
                if s_to > branch["rate_a"] && s_to - branch["rate_a"] > sm_vio
                    sm_vio = s_to - branch["rate_a"]
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


"returns a sorted list of branch flow violations"
function branch_violations_sorted_ratec(network, solution)
    branch_violations = []

    if haskey(solution, "branch")
        for (i,branch) in network["branch"]
            if branch["br_status"] != 0
                branch_sol = solution["branch"][i]
                s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)

                sm_vio = 0.0
                # TODO update to rate_c
                if s_fr > branch["rate_c"]
                    sm_vio = s_fr - branch["rate_c"]
                end
                if s_to > branch["rate_c"] && s_to - branch["rate_c"] > sm_vio
                    sm_vio = s_to - branch["rate_c"]
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



"assumes a vaild ac solution is included in the data and computes the branch flow values"
function calc_branch_flow_ac_goc(data::Dict{String,<:Any})
    @assert("per_unit" in keys(data) && data["per_unit"])
    @assert(!haskey(data, "conductors"))

    if InfrastructureModels.ismultinetwork(data)
        nws = Dict{String,Any}()
        for (i,nw_data) in data["nw"]
            nws[i] = _calc_branch_flow_ac_goc(nw_data)
        end
        return Dict{String,Any}(
            "nw" => nws,
            "per_unit" => data["per_unit"],
            "baseMVA" => data["baseMVA"]
        )
    else
        flows = _calc_branch_flow_ac_goc(data)
        flows["per_unit"] = data["per_unit"]
        flows["baseMVA"] = data["baseMVA"]
        return flows
    end
end


"helper function for calc_branch_flow_ac"
function _calc_branch_flow_ac_goc(data::Dict{String,<:Any})
    vm = Dict(bus["index"] => bus["vm"] for (i,bus) in data["bus"])
    va = Dict(bus["index"] => bus["va"] for (i,bus) in data["bus"])

    flows = Dict{String,Any}()
    for (i,branch) in data["branch"]
        if branch["br_status"] != 0
            f_bus = branch["f_bus"]
            t_bus = branch["t_bus"]

            g, b = calc_branch_y(branch)
            tr, ti = calc_branch_t(branch)
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

function compute_power_balance_deltas!(network)
    flows = calc_branch_flow_ac_goc(network)
    PowerModels.update_data!(network, flows)

    balance = PowerModels.calc_power_balance(network)
    PowerModels.update_data!(network, balance)

    p_delta_abs = [abs(bus["p_delta"]) for (i,bus) in network["bus"] if bus["bus_type"] != 4]
    q_delta_abs = [abs(bus["q_delta"]) for (i,bus) in network["bus"] if bus["bus_type"] != 4]

    return (
        p_delta_abs_max = maximum(p_delta_abs),
        p_delta_abs_mean = mean(p_delta_abs),
        q_delta_abs_max = maximum(q_delta_abs),
        q_delta_abs_mean = mean(q_delta_abs),
    )
end


##### GOC Solution Writers #####

function build_contingency_solutions(network, solution)
    contingency_solutions = Dict{String,Any}()

    solution = deepcopy(solution)
    solution["delta"] = 0.0

    for cont in network["gen_contingencies"]
        contingency_solutions[cont.label] = solution
    end

    for cont in network["branch_contingencies"]
        contingency_solutions[cont.label] = solution
    end

    return contingency_solutions
end




"checks feasibility criteria of network solution, produces an error if a problem is found"
function check_network_solution(network)
    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            if bus["vm"] > bus["vmax"] || bus["vm"] < bus["vmin"]
                error(LOGGER, "vm on $(bus["source_id"]) is not in bounds $(bus["vmin"]) to $(bus["vmax"]), given $(bus["vm"])")
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
                        error(LOGGER, "bs on $(shunt["source_id"]) is not in bounds $(shunt["bmin"]) to $(shunt["bmax"]), given $(shunt["bs"])")
                    end
                end
            end
        end
    end

    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            if gen["pg"] > gen["pmax"] || gen["pg"] < gen["pmin"]
                error(LOGGER, "pg on gen $(gen["source_id"]) not in bounds $(gen["pmin"]) to $(gen["pmax"]), given $(gen["pg"])")
            end

            if gen["qg"] > gen["qmax"] || gen["qg"] < gen["qmin"]
                error(LOGGER, "pg on gen $(gen["source_id"]) not in bounds $(gen["qmin"]) to $(gen["qmax"]), given $(gen["qg"])")
            end
        end
    end
end


"checks feasibility criteria of network solution, corrects when possible"
function correct_network_solution!(network)

    # default value is required for correctness of max and mean computations when no changes are made 
    vm_changes = [0.0]
    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            if bus["vm"] > bus["vmax"]
                warn(LOGGER, "update vm on bus $(i) to be in bounds $(bus["vm"]) -> $(bus["vmax"])")
                push!(vm_changes, bus["vm"] - bus["vmax"])
                bus["vm"] = bus["vmax"]
            end

            if bus["vm"] < bus["vmin"]
                warn(LOGGER, "update vm on bus $(i) to be in bounds $(bus["vm"]) -> $(bus["vmin"])")
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
                    warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmax"])")
                    push!(bs_changes, shunt["bs"] - shunt["bmax"])
                    shunt["bs"] = shunt["bmax"]
                end
                if shunt["bs"] < shunt["bmin"]
                    warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmin"])")
                    push!(bs_changes, shunt["bmin"] - shunt["bs"])
                    shunt["bs"] = shunt["bmin"]
                end
            end
        else
            warn(LOGGER, "shunt $(i) missing dispatchable parameter")
        end
    end

    pg_changes = [0.0]
    qg_changes = [0.0]
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            if gen["pg"] > gen["pmax"]
                warn(LOGGER, "update pg on gen $(i) to be in bounds $(gen["pg"]) -> $(gen["pmax"])")
                push!(pg_changes, gen["pg"] - gen["pmax"])
                gen["pg"] = gen["pmax"]
            end
            if gen["pg"] < gen["pmin"]
                warn(LOGGER, "update pg on gen $(i) to be in bounds $(gen["pg"]) -> $(gen["pmin"])")
                push!(pg_changes, gen["pmin"] - gen["pg"])
                gen["pg"] = gen["pmin"]
            end

            if gen["qg"] > gen["qmax"]
                warn(LOGGER, "update qg on gen $(i) to be in bounds $(gen["qg"]) -> $(gen["qmax"])")
                push!(qg_changes, gen["qg"] - gen["qmax"])
                gen["qg"] = gen["qmax"]
            end
            if gen["qg"] < gen["qmin"]
                warn(LOGGER, "update qg on gen $(i) to be in bounds $(gen["qg"]) -> $(gen["qmin"])")
                push!(qg_changes, gen["qmin"] - gen["qg"])
                gen["qg"] = gen["qmin"]
            end
        else
            gen["pg"] = 0.0
            gen["qg"] = 0.0
        end
    end

    _summary_changes(network, "base_case", vm_changes, bs_changes, pg_changes, qg_changes)
end



"checks feasibility criteria of contingencies, corrects when possible"
function correct_contingency_solutions!(network, contingency_solutions)
    bus_gens = gens_by_bus(network)

    cont_changes = Int64[]
    cont_vm_changes_max = [0.0]
    cont_bs_changes_max = [0.0]
    cont_pg_changes_max = [0.0]
    cont_qg_changes_max = [0.0]

    for cont_sol in contingency_solutions
        changes = correct_contingency_solution!(network, cont_sol; bus_gens=bus_gens)

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
function correct_contingency_solution!(network, cont_sol; bus_gens = gens_by_bus(network))
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
                    if qg >= qmax && bus["vm"] - vm_eq_tol/10 > nw_bus["vm"]
                        #println(qmin, " ", qg, " ", qmax)
                        #warn(LOGGER, "update qg $(qmin) $(qg) $(qmax)")
                        warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match set-point $(bus["vm"]) -> $(nw_bus["vm"]) due to qg upper bound and vm direction")
                        push!(vm_changes, abs(bus["vm"] - nw_bus["vm"]))
                        bus["vm"] = nw_bus["vm"]
                    end

                    if qg <= qmin && bus["vm"] + vm_eq_tol/10 < nw_bus["vm"]
                        #println(qmin, " ", qg, " ", qmax)
                        #warn(LOGGER, "update qg $(qmin) $(qg) $(qmax)")
                        warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match set-point $(bus["vm"]) -> $(nw_bus["vm"]) due to qg lower bound and vm direction")
                        push!(vm_changes, abs(bus["vm"] - nw_bus["vm"]))
                        bus["vm"] = nw_bus["vm"]
                    end
                end
            end

            if bus["vm"] > nw_bus["vmax"]
                warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match ub $(bus["vm"]) -> $(nw_bus["vmax"]) due to out of bounds")
                push!(vm_changes, bus["vm"] - nw_bus["vmax"])
                bus["vm"] = nw_bus["vmax"]
            end

            if bus["vm"] < nw_bus["vmin"]
                warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match lb $(bus["vm"]) -> $(nw_bus["vmin"]) due to out of bounds")
                push!(vm_changes, nw_bus["vmin"] - bus["vm"])
                bus["vm"] = nw_bus["vmin"]
            end
        else
            bus["vm"] = 0.0
            bus["va"] = 0.0
        end
    end


    bs_changes = [0.0]
    for (i,shunt) in cont_sol["shunt"]
        nw_shunt = network["shunt"][i]
        if haskey(nw_shunt, "dispatchable") && nw_shunt["dispatchable"]
            @assert nw_shunt["gs"] == 0.0
            @assert haskey(nw_shunt, "bmin") && haskey(nw_shunt, "bmax")
            if shunt["bs"] > nw_shunt["bmax"]
                warn(LOGGER, "update bs on shunt $(i) in contingency $(label) to be in bounds $(shunt["bs"]) -> $(nw_shunt["bmax"])")
                push!(bs_changes, shunt["bs"] - nw_shunt["bmax"])
                shunt["bs"] = nw_shunt["bmax"]
            end
            if shunt["bs"] < nw_shunt["bmin"]
                warn(LOGGER, "update bs on shunt $(i) in contingency $(label) to be in bounds $(shunt["bs"]) -> $(nw_shunt["bmin"])")
                push!(bs_changes, nw_shunt["bmin"] - shunt["bs"])
                shunt["bs"] = nw_shunt["bmin"]
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
                if !isapprox(bus["vm"], nw_bus["vm"], atol=vm_eq_tol/2)
                    #debug(LOGGER, "bus $(bus_id) : vm_base $(nw_bus["vm"]) - vm $(bus["vm"]) : reactive bounds $(nw_gen["qmin"]) - $(gen["qg"]) - $(nw_gen["qmax"])")
                    warn(LOGGER, "update vm on bus $(bus_id) in contingency $(label) to match base case $(bus["vm"]) -> $(nw_bus["vm"]) due to within reactive bounds")
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
                warn(LOGGER, "pg value on gen $(i) $(nw_gen["source_id"]) in contingency $(label) is not consistent with the computed value given:$(gen["pg"]) calc:$(pg_calc)")
            end

            if gen["pg"] > nw_gen["pmax"]
                warn(LOGGER, "update pg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match ub $(gen["pg"]) -> $(nw_gen["pmax"])")
                push!(pg_changes, gen["pg"] - nw_gen["pmax"])
                gen["pg"] = nw_gen["pmax"]
            end

            if gen["pg"] < nw_gen["pmin"]
                warn(LOGGER, "update pg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match lb $(gen["pg"]) -> $(nw_gen["pmin"])")
                push!(pg_changes, nw_gen["pmin"] - gen["pg"])
                gen["pg"] = nw_gen["pmin"]
            end

            if gen["qg"] > nw_gen["qmax"]
                warn(LOGGER, "update qg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match ub $(gen["qg"]) -> $(nw_gen["qmax"])")
                push!(qg_changes, gen["qg"] - nw_gen["qmax"])
                gen["qg"] = nw_gen["qmax"]
            end

            if gen["qg"] < nw_gen["qmin"]
                warn(LOGGER, "update qg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match lb $(gen["qg"]) -> $(nw_gen["qmin"])")
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
    PowerModels.update_data!(network_tmp, cont_sol)
    deltas = compute_power_balance_deltas!(network_tmp)
    println("$(label): $(deltas)")
    =#

    cont_changed = length(vm_changes) > 1 || length(bs_changes) > 1 || length(pg_changes) > 1 || length(qg_changes) > 1

    if cont_changed
        _summary_changes(network, label, vm_changes, bs_changes, pg_changes, qg_changes)
    end

    return (changed=Int(cont_changed), vm_changes_max=maximum(vm_changes), bs_changes_max=maximum(bs_changes), pg_changes_max=maximum(pg_changes), qg_changes_max=maximum(qg_changes))
end


function _summary_changes(network, contingency, vm_changes, bs_changes, pg_changes, qg_changes)
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
    info(LOGGER, join(data, ", "))

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
    info(LOGGER, join(data, ", "))
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




function write_solution(goc_data, pm_network, contingencies; solution_1="", solution_2="")
    files = goc_data.files
    if length(solution_1) > 0
        files["sol1"] = solution_1
    else
        ini_dir = dirname(goc_data.ini_file)
        files["sol1"] = joinpath(ini_dir, goc_data.scenario, "solution1.txt")
    end

    if length(solution_2) > 0
        files["sol2"] = solution_2
    else
        ini_dir = dirname(goc_data.ini_file)
        files["sol2"] = joinpath(ini_dir, goc_data.scenario, "solution2.txt")
    end

    open("files.json", "w") do input
        JSON.print(input, files)
    end

    write_solution1(pm_network, solution_file=files["sol1"])
    write_solution2(pm_network, contingencies, solution_file=files["sol2"])
end


function write_solution1(pm_network; output_dir="", solution_file="solution1.txt")
    if length(output_dir) > 0
        solution_path = joinpath(output_dir, solution_file)
    else
        solution_path = solution_file
    end

    open(solution_path, "w") do sol1
        base_mva = pm_network["baseMVA"]

        bus_switched_shunt_b = Dict(i => 0.0 for (i,bus) in pm_network["bus"])
        for (i,shunt) in pm_network["shunt"]
            # test checks if a "switched shunt" in the orginal data model
            if shunt["dispatchable"] && shunt["status"] == 1
                @assert shunt["gs"] == 0.0
                bus_switched_shunt_b["$(shunt["shunt_bus"])"] += shunt["bs"]
            end
        end

        write(sol1, "-- bus section\n")
        write(sol1, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
        for (i,bus) in pm_network["bus"]
            write(sol1, "$(bus["index"]), $(bus["vm"]), $(rad2deg(bus["va"])), $(base_mva*bus_switched_shunt_b[i])\n")
        end

        write(sol1, "-- generator section\n")
        write(sol1, "i, id, p(MW), q(MVAR)\n")
        for (i,gen) in pm_network["gen"]
            bus_index = gen["source_id"][2]
            gen_id = gen["source_id"][3]
            write(sol1, "$(bus_index), $(gen_id), $(base_mva*gen["pg"]), $(base_mva*gen["qg"])\n")
        end
    end
end


function write_solution2(pm_network, contingencies; output_dir="", solution_file="solution2.txt")
    if length(output_dir) > 0
        solution_path = joinpath(output_dir, solution_file)
    else
        solution_path = solution_file
    end

    open(solution_path, "w") do sol2
        base_mva = pm_network["baseMVA"]

        for cont_solution in contingencies
            write_solution2_contingency(sol2, pm_network, cont_solution)
        end
    end

    return solution_path
end


function write_solution2_contingency(io::IO, pm_network, contingency_solution)
    base_mva = pm_network["baseMVA"]

    bus_switched_shunt_b = Dict(i => 0.0 for (i,bus) in pm_network["bus"])
    for (i,nw_shunt) in pm_network["shunt"]
        if nw_shunt["dispatchable"] && nw_shunt["status"] == 1
            #@assert nw_shunt["gs"] == 0.0
            shunt = contingency_solution["shunt"][i]
            bus_switched_shunt_b["$(nw_shunt["shunt_bus"])"] += shunt["bs"]
        end
    end

    write(io, "-- contingency\n")
    write(io, "label\n")
    write(io, "$(contingency_solution["label"])\n")

    write(io, "-- bus section\n")
    write(io, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
    for (i,bus) in contingency_solution["bus"]
        nw_bus = pm_network["bus"][i]
        write(io, "$(nw_bus["index"]), $(bus["vm"]), $(rad2deg(bus["va"])), $(base_mva*bus_switched_shunt_b[i])\n")
    end

    write(io, "-- generator section\n")
    write(io, "i, id, p(MW), q(MVAR)\n")
    for (i,gen) in contingency_solution["gen"]
        nw_gen = pm_network["gen"][i]
        bus_index = nw_gen["source_id"][2]
        gen_id = nw_gen["source_id"][3]
        write(io, "$(bus_index), $(gen_id), $(base_mva*gen["pg"]), $(base_mva*gen["qg"])\n")
    end

    write(io, "-- delta section\n")
    write(io, "delta(MW)\n")
    write(io, "$(base_mva*contingency_solution["delta"])\n")
end


function combine_files(files, output_file_name; output_dir="")
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


function remove_files(files)
    for file in files
        if isfile(file)
            info(LOGGER, "deleting: $(file)")
            rm(file)
        else
            info(LOGGER, "skipping file: $(file)")
        end
    end
end


function remove_solution_files(;output_dir="", solution1_file="solution1.txt", solution2_file="solution2.txt")
    if length(output_dir) > 0
        solution1_path = joinpath(output_dir, solution1_file)
    else
        solution1_path = solution1_file
    end

    if isfile(solution1_path)
        info(LOGGER, "deleting: $(solution1_path)")
        rm(solution1_path)
    else
        info(LOGGER, "skipping file: $(solution1_path)")
    end


    if length(output_dir) > 0
        solution2_path = joinpath(output_dir, solution2_file)
    else
        solution2_path = solution2_file
    end

    if isfile(solution2_path)
        info(LOGGER, "deleting: $(solution2_path)")
        rm(solution2_path)
    else
        info(LOGGER, "skipping file: $(solution2_path)")
    end
end

function remove_detail_file(;output_dir="", detail_file="detail.csv")
    if length(output_dir) > 0
        detail_file_path = joinpath(output_dir, detail_file)
    else
        detail_file_path = detail_file
    end

    if isfile(detail_file_path)
        info(LOGGER, "deleting: $(detail_file_path)")
        rm(detail_file_path)
    else
        info(LOGGER, "skipping file: $(detail_file_path)")
    end
end


function write_file_paths(files; output_dir="", solution1_file="solution1.txt", solution2_file="solution2.txt")
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


function write_scopf_summary(scenario_id, pm_network, objective; branch_flow_cuts=0, objective_lb=-Inf, load_time=-1.0, solve_time=-1.0, filter_time=-1.0, total_time=load_time+solve_time+filter_time)
    println("")

    data = [
        "----",
        "scenario id",
        "num bus",
        "num branch",
        "num gen",
        "num load",
        "num shunt",
        "num gen cont",
        "num branch cont",
        "active gen cont",
        "active branch cont",
        "branch flow cuts",
        "objective ub",
        "objective lb",
        "load time (sec.)",
        "solve time (sec.)",
        "filter time (sec.)",
        "total time (sec.)",
    ]
    println(join(data, ", "))

    data = [
        "DATA",
        scenario_id,
        length(pm_network["bus"]),
        length(pm_network["branch"]),
        length(pm_network["gen"]),
        length(pm_network["load"]),
        length(pm_network["shunt"]),
        length(pm_network["gen_contingencies"]),
        length(pm_network["branch_contingencies"]),
        length(pm_network["gen_contingencies_active"]),
        length(pm_network["branch_contingencies_active"]),
        branch_flow_cuts,
        objective,
        objective_lb,
        load_time,
        solve_time,
        filter_time,
        total_time,
    ]
    println(join(data, ", "))

end


function write_power_balance_summary(scenario_id, p_delta_abs_max, q_delta_abs_max, p_delta_abs_mean, q_delta_abs_mean)
    println("")

    data = [
        "----",
        "scenario id",
        "p_delta abs max",
        "q_delta abs max",
        "p_delta abs mean",
        "q_delta abs mean",
    ]
    println(join(data, ", "))

    data = [
        "DATA_PB",
        scenario_id,
        p_delta_abs_max,
        q_delta_abs_max,
        p_delta_abs_mean,
        q_delta_abs_mean,
    ]
    println(join(data, ", "))

end
