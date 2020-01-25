

@testset "pg response" for case in cases

    network_s1 = build_pm_model(case)
    solution1_file = joinpath(dirname(case.files["raw"]), "solution1.txt")
    pm_sol = build_pm_solution(network_s1, solution1_file)
    PowerModels.update_data!(network_s1, pm_sol)


    for (i,gen) in network_s1["gen"]
        gen["pg_base"] = gen["pg"]
    end

    pg_target = 0.0
    for (i,gen) in network_s1["gen"]
        if gen["gen_status"] != 0
            pg_target += gen["pg"]
        end
    end


    for (i,gen) in network_s1["gen"]
        if gen["gen_status"] == 0
            continue
        end

        gen_bus = network_s1["bus"]["$(gen["gen_bus"])"]
        area_gens = network_s1["area_gens"][gen_bus["area"]]


        # NOTE this is a slow operation
        network_cont = deepcopy(network_s1)

        cont_gen = network_cont["gen"][i]
        cont_gen["gen_status"] = 0
        pg_missing = cont_gen["pg"]

        #println("")
        #info(LOGGER, "generation removed: $(pg_missing)")

        apply_pg_response!(network_cont, area_gens, pg_missing)

        pg_total = 0.0
        for (i,gen) in network_cont["gen"]
            if gen["gen_status"] != 0
                pg_total += gen["pg"]
            end
        end

        pg_comp = comp_pg_response_total(network_cont, area_gens)

        #info(LOGGER, "delta $delta")
        #info(LOGGER, "gen $(pg_target) / $(pg_total) / $(pg_comp)")

        @test isapprox(pg_target, pg_total)
        @test isapprox(pg_total, pg_comp)
    end

end

