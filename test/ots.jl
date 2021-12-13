# Tests formulations in prob/ots
@testset "test ots" begin

@testset "c2 ots" begin

    opf_ac_objective = [593064.85, 379600.54, 543101.32]
    @testset "ots acp - $(i)" for (i,network) in enumerate(c2_networks)
        #parse_c2_opf_files(c1_ini_file, scenario_id=c1_scenarios[i])
        result = run_c2_ots_soft_bus(network, ACPPowerModel, nlp_solver, relax_integrality=true)

        @test isapprox(result["termination_status"], LOCALLY_SOLVED)
        @test isapprox(result["objective"], opf_ac_objective[i]; atol = 1e0)
    end

end

end # close test group
