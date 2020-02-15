# Tests formulations in prob/opf

opf_cheap_dc_objective = [14642.16, 26982.17]
@testset "opf cheap dc - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_cheap_dc(network, DCPPowerModel, lp_solver)

    @test isapprox(result["termination_status"], OPTIMAL)
    @test isapprox(result["objective"], opf_cheap_dc_objective[i]; atol = 1e0)
end


opf_pg_pf_rect_objective = [123117.35405092733, 9.480588301784407e6]
@testset "opf pg pf rect 5 - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)

    deactivate_rate_a!(network)
    activate_rate_a_violations!(network)

    result = run_opf_pg_pf_rect_5(network, ACRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_pg_pf_rect_objective[i]; atol = 1e0)
end


