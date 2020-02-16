# Tests formulations in prob/security-stage

function security_stage_data(network)
    network = deepcopy(network)

    network["cont_label"] = "testing"
    network["delta"] = 0.0

    for (a,gens) in network["area_gens"]
        network["response_gens"] = gens
        break
    end

    for (i,branch) in network["branch"]
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        branch["g"] = g
        branch["b"] = b
        branch["tr"] = tr
        branch["ti"] = ti
    end

    bus_gens = gens_by_bus(network)

    network["delta"] = 0.0
    for (i,bus) in network["bus"]
        bus["vm_base"] = bus["vm"]
        bus["vm_start"] = bus["vm"]
        bus["va_start"] = bus["va"]
        bus["vm_fixed"] = length(bus_gens[i]) != 0
    end

    for (i,gen) in network["gen"]
        gen["pg_base"] = gen["pg"]
        gen["pg_start"] = gen["pg"]
        gen["qg_start"] = gen["qg"]
        gen["pg_fixed"] = false
        gen["qg_fixed"] = false
    end

    return network
end


@testset "pf soft rect - $(i)" for (i,network) in enumerate(networks)
    network = security_stage_data(network)

    result = run_pf_soft_rect(network, ACRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf soft rect fixedpoint - $(i)" for (i,network) in enumerate(networks)
    network = security_stage_data(network)
    pg_lost = 0.0

    result = run_fixpoint_pf_soft!(network, pg_lost, ACRPowerModel, nlp_solver, iteration_limit=10)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf soft rect fixedpoint - infeasible" begin
    network = security_stage_data(network_infeasible)
    pg_lost = 0.0

    result = run_fixpoint_pf_soft!(network, pg_lost, ACRPowerModel, nlp_solver, iteration_limit=10)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end


@testset "pf fixed nbf rect2 - $(i)" for (i,network) in enumerate(networks)
    network = security_stage_data(network)

    result = run_fixed_pf_nbf_rect2(network, ACRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf fixed nbf rect2 ds- $(i)" for (i,network) in enumerate(networks)
    network = security_stage_data(network)

    for (i,gen) in network["gen"]
        if gen["index"] in network["response_gens"]
            gen["pg_fixed"] = true
        end
    end

    result = run_fixed_pf_nbf_rect2_ds(network, ACRPowerModel, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf v2_3 fixpoint - $(i)" for (i,network) in enumerate(networks)
    network = security_stage_data(network)
    pg_lost = 0.0

    result = run_fixpoint_pf_v2_3!(network, pg_lost, ACRPowerModel, nlp_solver, iteration_limit=5)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf v2_3 fixpoint - infeasible" begin
    network = security_stage_data(network_infeasible)
    pg_lost = 0.0

    result = run_fixpoint_pf_v2_3!(network, pg_lost, ACRPowerModel, nlp_solver, iteration_limit=5)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
end
