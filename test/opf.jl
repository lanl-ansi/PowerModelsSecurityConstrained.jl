# Tests formulations in prob/opf
@testset "test opf" begin

opf_shunt_ac_objective = [14676.95, 27545.01]
@testset "opf shunt acp - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_shunt(network, ACPPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_shunt_ac_objective[i]; atol = 1e0)
end

@testset "opf shunt acr - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_shunt(network, ACRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_shunt_ac_objective[i]; atol = 1e0)
end

opf_shunt_soc_objective = [14666.91, 27510.25]
@testset "opf shunt soc - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_shunt(network, SOCWRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_shunt_soc_objective[i]; atol = 1e0)
end

opf_shunt_dcp_objective = [14642.16, 26982.17]
@testset "opf shunt dcp - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_shunt(network, DCPPowerModel, lp_solver)

    @test isapprox(result["termination_status"], OPTIMAL)
    @test isapprox(result["objective"], opf_shunt_dcp_objective[i]; atol = 1e0)
end

@testset "opf shunt dc - infeasible" begin
    result = run_opf_shunt(network_infeasible, DCPPowerModel, lp_solver)
    @test isapprox(result["termination_status"], INFEASIBLE)
end


opf_cheap_ac_objective = [14676.95, 27542.05]
@testset "opf cheap acp - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_cheap(network, ACPPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_cheap_ac_objective[i]; atol = 1e0)
end

@testset "opf cheap acr - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_cheap(network, ACRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_cheap_ac_objective[i]; atol = 1e0)
end

opf_cheap_soc_objective = [14666.81, 27507.30]
@testset "opf cheap soc - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_cheap(network, SOCWRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], opf_cheap_soc_objective[i]; atol = 1e0)
end

# not fully general
# opf_cheap_dcp_objective = [14676.95, 27542.05]
# @testset "opf cheap dcp - $(i)" for (i,network) in enumerate(networks)

#     result = run_opf_cheap(network, ACPPowerModel, lp_solver)

#     @test isapprox(result["termination_status"], OPTIMAL)
#     @test isapprox(result["objective"], opf_cheap_ac_objective[i]; atol = 1e0)
# end



opf_cheap_dc_objective = [14642.16, 26982.17]
@testset "opf cheap_dc dc - $(i)" for (i,network) in enumerate(networks)

    result = run_opf_cheap_dc(network, DCPPowerModel, lp_solver)

    @test isapprox(result["termination_status"], OPTIMAL)
    @test isapprox(result["objective"], opf_cheap_dc_objective[i]; atol = 1e0)
end

@testset "opf shunt dc - infeasible" begin
    result = run_opf_cheap_dc(network_infeasible, DCPPowerModel, lp_solver)
    @test isapprox(result["termination_status"], INFEASIBLE)
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

@testset "opf pg pf rect 5 - infeasible" begin
    network = deepcopy(network_infeasible)

    deactivate_rate_a!(network)
    activate_rate_a_violations!(network)

    result = run_opf_pg_pf_rect_5(network, ACRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_INFEASIBLE)
end

end # close test group