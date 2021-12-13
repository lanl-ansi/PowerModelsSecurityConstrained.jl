@testset "test writing and reading files" begin

@testset "c1 io" begin

    @testset "basic read/write process" begin
        pm_data = deepcopy(c1_networks[1])

        correct_c1_solution!(pm_data)
        write_c1_solution1(pm_data, solution_file="solution1-tmp.txt")

        pm_sol = read_c1_solution1(pm_data, state_file="solution1-tmp.txt")
        remove_c1_solution_files(solution1_file="solution1-tmp.txt")

        PowerModels.update_data!(pm_data, pm_sol)

        result = run_c1_opf_shunt(pm_data, ACPPowerModel, nlp_solver)

        @test isapprox(result["termination_status"], LOCALLY_SOLVED)
        @test isapprox(result["objective"], 14676.89; atol = 1e0)

        result["solution"]["label"] = "tmp"
        result["solution"]["cont_type"] = "branch"
        result["solution"]["cont_comp_id"] = 1
        result["solution"]["delta"] = 0.0

        correct_c1_contingency_solution!(pm_data, result["solution"])

        result["solution"]["label"] = "tmp"
        result["solution"]["cont_type"] = "gen"
        result["solution"]["cont_comp_id"] = 1
        result["solution"]["delta"] = 0.0

        correct_c1_contingency_solution!(pm_data, result["solution"])
    end

end


@testset "c2 io" begin

    @testset "basic read/write process" begin
        pm_data = deepcopy(c2_networks[1])

        normalize_tm_values!(pm_data)
        set_va_start_values!(pm_data)

        correct_c2_solution!(pm_data)
        write_c2_solution(pm_data, label="TMP")


        pm_sol = read_c2_solution(pm_data, label="TMP")
        remove_c2_solution_files()

        PowerModels.update_data!(pm_data, pm_sol)
        update_previous!(pm_data)

        set_va_start_values!(pm_data)

        result = run_c2_opf_soft_ctg(pm_data, ACPPowerModel, nlp_solver)

        @test isapprox(result["termination_status"], LOCALLY_SOLVED)
        @test isapprox(result["objective"], 593064.46; atol = 1e0)
    end

end

end