# Tests formulations in prob/scopf

scopf_dc_cuts_soft_woc_objective = [14642.16, 26982.17]
@testset "scopf cuts dc soft 2, without cuts - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)
    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = run_scopf_cuts_dc_soft_2(network, DCPPowerModel, lp_solver)

    @test isapprox(result["termination_status"], OPTIMAL)
    @test isapprox(result["objective"], scopf_dc_cuts_soft_woc_objective[i]; atol = 1e0)
end

@testset "scopf cuts dc soft 2, infeasible" begin
    network = deepcopy(network_infeasible)
    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = run_scopf_cuts_dc_soft_2(network, DCPPowerModel, lp_solver)
    @test isapprox(result["termination_status"], INFEASIBLE)
end

scopf_dc_cuts_soft_wc_objective = [14642.16, 30737.94]
@testset "scopf cuts dc soft 2, with cuts - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)
    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = run_opf_cheap_dc(network, DCPPowerModel, lp_solver)
    @test isapprox(result["termination_status"], OPTIMAL)

    update_active_power_data!(network, result["solution"])

    cuts = check_contingencies_branch_flow_ratec_nd_first_lazy(network, total_cut_limit=2, gen_flow_cuts=[], branch_flow_cuts=[])

    #println(length(cuts.gen_cuts) + length(cuts.branch_cuts))
    #cuts_found = sum(length(c.gen_cuts)+length(c.branch_cuts) for c in cuts)
    append!(network["gen_flow_cuts"], cuts.gen_cuts)
    append!(network["branch_flow_cuts"], cuts.branch_cuts)

    result = run_scopf_cuts_dc_soft_2(network, DCPPowerModel, lp_solver)

    @test isapprox(result["termination_status"], OPTIMAL)
    @test isapprox(result["objective"], scopf_dc_cuts_soft_wc_objective[i]; atol = 1e0)
end

# opf_pg_pf_rect_objective = [123117.35405092733, 9.480588301784407e6]
# @testset "opf pg pf rect 5" for (i,network) in enumerate(networks)
#     network = deepcopy(network)

#     deactivate_rate_a!(network)
#     activate_rate_a_violations!(network)

#     result = run_opf_pg_pf_rect_5(network, ACRPowerModel, nlp_solver)

#     @test isapprox(result["termination_status"], LOCALLY_SOLVED)
#     @test isapprox(result["objective"], opf_pg_pf_rect_objective[i]; atol = 1e0)
# end


