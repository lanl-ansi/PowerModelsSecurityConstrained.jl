using ArgParse
using PowerModelsSecurityConstrained

function scopf_main(parsed_args)
    ini_file = parsed_args["file"]
    scenario = parsed_args["scenario"]
    scoring_method = parsed_args["scoring-method"]
    time_limit_sol1 = parsed_args["time-limit"]
    time_limit_total = parsed_args["hard-time-limit"]
    distribute = parsed_args["distribute"]
    skip_solution2 = parsed_args["skip-solution2"]
    remove_solutions = parsed_args["remove-solutions"]
    gurobi = parsed_args["gurobi"]

    files, scenario_id = find_goc_files(ini_file, scenario_id=scenario)

    ini_dir = dirname(ini_file)
    output_dir = joinpath(ini_dir, scenario_id)

    info(LOGGER, "output dir: $(output_dir)")

    time_sol1_start = time()
    compute_solution1(files["con"], files["inl"], files["raw"], files["rop"], time_limit_sol1-60, scoring_method, "network name"; output_dir=output_dir, scenario_id=scenario_id, gurobi=gurobi)
    time_sol1 = time() - time_sol1_start

    if !skip_solution2
        compute_solution2(files["con"], files["inl"], files["raw"], files["rop"], trunc(Int, time_limit_total-time_sol1-60), scoring_method, "network name"; output_dir=output_dir, scenario_id=scenario_id)
    end

    if remove_solutions
        remove_solution_files(output_dir=output_dir)
        remove_detail_file(output_dir=output_dir)
    end
end


function parse_scopf_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file", "-f"
            help = "the data initiation file (.ini)"
            required = true
        "--scenario", "-s"
            help = "the scenario to run (directory)"
            default = ""
        "--scoring-method", "-m"
            help = "the objective function (1,2,3,4)"
            arg_type = Int
            default = 2
        "--time-limit", "-t"
            help = "solution 1 runtime limit (seconds)"
            arg_type = Int
            default = 600000 #2700
        "--hard-time-limit"
            help = "solution 1 and 2 runtime limit (seconds)"
            arg_type = Int
            default = 600000 #34200
        "--distribute", "-d"
            help = "run on multiple processes"
            action = :store_true
        "--skip-solution2"
            help = "skip computation of solution2 file"
            action = :store_true
        "--remove-solutions"
            help = "delete the solution files after competition"
            action = :store_true
        "--gurobi", "-g"
            help = "use Gurobi for solving lp, qp and mip problems"
            action = :store_true
            #default = false
    end

    return parse_args(s)
end
