# Tests formulations in prob/security-stage
@testset "test security-stage" begin

function contingency_stage_data(network)
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


solution2_lines = [31,600]
@testset "write_solution2 - $(i)" for (i,network) in enumerate(networks)
    network = contingency_stage_data(network)

    bus_gens = gens_by_bus(network)
    cont = network["gen_contingencies"][1]

    cont_gen = network["gen"]["$(cont.idx)"]
    cont_gen["contingency"] = true
    cont_gen["gen_status"] = 0
    pg_lost = cont_gen["pg"]

    gen_bus = network["bus"]["$(cont_gen["gen_bus"])"]
    if length(bus_gens["$(gen_bus["index"])"]) == 1
        gen_bus["vm_fixed"] = false
    end

    result = run_fixpoint_pf_bqv!(network, pg_lost, nlp_solver)
    #result = run_fixpoint_pf_bqv_native!(network, pg_lost, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)

    cont_sol = result["solution"]
    cont_sol["label"] = cont.label
    cont_sol["feasible"] = (result["termination_status"] == LOCALLY_SOLVED)
    cont_sol["cont_type"] = "gen"
    cont_sol["cont_comp_id"] = cont.idx

    cont_sol["gen"]["$(cont.idx)"] = Dict("pg" => 0.0, "qg" => 0.0)
    cont_sol["delta"] = 0.0

    write_solution2(network, [cont_sol], solution_file="solution2-tmp.txt")

    open("solution2-tmp.txt", "r") do file
        lines = countlines(file)
        @test lines == solution2_lines[i]
    end
    rm("solution2-tmp.txt")
end


@testset "pf soft rect - $(i)" for (i,network) in enumerate(networks)
    network = contingency_stage_data(network)

    result = run_pf_bqv_acr(network, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf soft rect fixedpoint - $(i)" for (i,network) in enumerate(networks)
    network = contingency_stage_data(network)
    pg_lost = 0.0

    result = run_fixpoint_pf_bqv!(network, pg_lost, nlp_solver, iteration_limit=10)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf soft rect fixedpoint - infeasible" begin
    network = contingency_stage_data(network_infeasible)
    pg_lost = 0.0

    result = run_fixpoint_pf_bqv!(network, pg_lost, nlp_solver, iteration_limit=10)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end


@testset "pf fixed nbf rect2 - $(i)" for (i,network) in enumerate(networks)
    network = contingency_stage_data(network)

    result = run_pf_fixed_acr(network, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf fixed nbf rect2 ds- $(i)" for (i,network) in enumerate(networks)
    network = contingency_stage_data(network)

    for (i,gen) in network["gen"]
        if gen["index"] in network["response_gens"]
            gen["pg_fixed"] = true
        end
    end

    result = run_pf_fixed_bp_slack_acr(network, nlp_solver)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf v2_3 fixpoint - $(i)" for (i,network) in enumerate(networks)
    network = contingency_stage_data(network)
    pg_lost = 0.0

    result = run_fixpoint_pf_pvpq!(network, pg_lost, nlp_solver, iteration_limit=5)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
    @test isapprox(result["objective"], 0.0; atol = 1e0)
end

@testset "pf v2_3 fixpoint - infeasible" begin
    network = contingency_stage_data(network_infeasible)
    pg_lost = 0.0

    result = run_fixpoint_pf_pvpq!(network, pg_lost, nlp_solver, iteration_limit=5)

    @test isapprox(result["termination_status"], LOCALLY_SOLVED)
end

end # close test group