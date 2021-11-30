@testset "test contigency" begin

@testset "pg response - $(i)" for (i,network) in enumerate(networks)

    network = deepcopy(network)
    PowerModels.update_data!(network, solutions[i])

    for (i,gen) in network["gen"]
        gen["pg_base"] = gen["pg"]
    end

    pg_target = 0.0
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            pg_target += gen["pg"]
        end
    end


    for (i,gen) in network["gen"]
        if gen["gen_status"] == 0
            continue
        end

        gen_bus = network["bus"]["$(gen["gen_bus"])"]
        area_gens = network["area_gens"][gen_bus["area"]]


        # NOTE this is a slow operation
        network_cont = deepcopy(network)

        cont_gen = network_cont["gen"][i]
        cont_gen["gen_status"] = 0
        pg_missing = cont_gen["pg"]

        #println("")
        #println("generation removed: $(pg_missing)")

        apply_c1_pg_response!(network_cont, area_gens, pg_missing)

        pg_total = 0.0
        for (i,gen) in network_cont["gen"]
            if gen["gen_status"] != 0
                pg_total += gen["pg"]
            end
        end

        pg_comp = calc_c1_pg_response_total(network_cont, area_gens)

        #println("delta $delta")
        #println("gen $(pg_target) / $(pg_total) / $(pg_comp)")

        @test isapprox(pg_target, pg_total)
        @test isapprox(pg_total, pg_comp)
    end

end



cuts_ratec_nd_first_lazy_gen = [0, 10]
cuts_ratec_nd_first_lazy_branch = [0, 10]
@testset "cuts ratec_nd_first_lazy - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)
    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = run_c1_opf_cheap(network, DCPPowerModel, lp_solver)
    @test isapprox(result["termination_status"], OPTIMAL)
    update_active_power_data!(network, result["solution"])

    cuts = check_c1_contingencies_branch_power_bpv(network, total_cut_limit=1000, gen_flow_cuts=[], branch_flow_cuts=[])

    @test isapprox(length(cuts.gen_cuts), cuts_ratec_nd_first_lazy_gen[i])
    @test isapprox(length(cuts.branch_cuts), cuts_ratec_nd_first_lazy_branch[i])
end


cuts_ratec_gen = [0, 10]
cuts_ratec_branch = [0, 10]
@testset "cuts ratec - $(i)" for (i,network) in enumerate(networks)
    network = deepcopy(network)
    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = run_c1_opf_cheap(network, DCPPowerModel, lp_solver)
    @test isapprox(result["termination_status"], OPTIMAL)
    update_active_power_data!(network, result["solution"])

    cuts = check_c1_contingencies_branch_power(network, total_cut_limit=1000, gen_flow_cuts=[], branch_flow_cuts=[])

    @test isapprox(length(cuts.gen_cuts), cuts_ratec_gen[i])
    @test isapprox(length(cuts.branch_cuts), cuts_ratec_branch[i])
end

end # close test group

