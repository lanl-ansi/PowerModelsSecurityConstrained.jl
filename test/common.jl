networks_opf = [build_pm_opf_model(case) for case in cases]
opf_ac_objective = [14676.9, 27564.91]
@testset "opf ac - $(i)" for (i,network) in enumerate(networks_opf)

    result = run_ac_opf(network, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_ac_objective[i]; atol = 1e0)
end


@testset "set_start_values - $(i)" for (i,network) in enumerate(networks_opf)
    network = deepcopy(network)

    set_start_values!(network)
    result = run_ac_opf(network, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_ac_objective[i]; atol = 1e0)
end


opf_ac_tight_objective = [14676.9, 28347.12]
@testset "opf ac tight- $(i)" for (i,network) in enumerate(networks_opf)
    network = deepcopy(network)

    tighten_constraints!(network)
    result = run_ac_opf(network, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_ac_tight_objective[i]; atol = 1e0)
end


mn_size = [2, 20]
@testset "build multinetwork - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)

    network["gen_contingencies"] = network["gen_contingencies"][1:min(10, length(network["gen_contingencies"]))]
    network["branch_contingencies"] = network["branch_contingencies"][1:min(10, length(network["branch_contingencies"]))]

    mn_case = build_scopf_multinetwork(network)

    @test isapprox(length(mn_case["nw"]), mn_size[i]; atol = 1e0)
end


@testset "rate_a deactivate / activate - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)

    deactivate_rate_a!(network)
    @test all(!haskey(branch, "rate_a") for (i,branch) in network["branch"])

    activate_rate_a!(network)
    @test all(haskey(branch, "rate_a") for (i,branch) in network["branch"])
end

first_cont_id = [9, 203]
@testset "contingency_order - $(i)" for (i,network) in enumerate(networks)
    order = contingency_order(network)
    @test isapprox(order[i].idx, first_cont_id[i]; atol = 1e0)
end

ref_bus_id = ["1", "17"]
@testset "correct_voltage_angles - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)
    update_data!(network, solutions[i])

    correct_voltage_angles!(network)

    @test isapprox(network["bus"][ref_bus_id[i]]["va"], 0.0; atol = 1e-2)
end


opf_p_delta_abs_max = [0.26746780838927353, 0.07392932526496376]
opf_q_delta_abs_max = [0.41571449226280965, 1.1036712386853447]
solution1_lines = [23,594]
@testset "write_solution1 - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)

    #deactivate_rate_a!(network)

    result = run_dc_opf(network, lp_solver)
    #result = run_dc_opf(network, lp_solver)
    @test isapprox(result["termination_status"], OPTIMAL)

    update_active_power_data!(network, result["solution"])

    balance = compute_power_balance_deltas!(network)
    @test isapprox(balance.p_delta_abs_max, opf_p_delta_abs_max[i])
    @test isapprox(balance.q_delta_abs_max, opf_q_delta_abs_max[i])

    correct_network_solution!(network)
    write_solution1(network, solution_file="solution1-tmp.txt")

    open("solution1-tmp.txt", "r") do file
        lines = countlines(file)
        @test lines == solution1_lines[i]
    end
    rm("solution1-tmp.txt")
end





