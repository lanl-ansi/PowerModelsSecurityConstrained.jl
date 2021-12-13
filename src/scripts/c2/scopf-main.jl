
using ArgParse

function scopf_c2_main(args)
    files = goc_c2_parse_case_files(args["case"])

    scoring_method = 4

    println("case size: $(countlines(files.raw))")

    println("running code1 and code2 with:")
    println("  $(files.raw)")
    println("  $(files.con)")
    println("  $(files.json)")
    println("  $("")")
    println("  $(args["time-limit"])")
    println("  $(scoring_method)")
    println("  $(files.case_id)")
    println("  $(files.scenario_id)")

    time_start = time()
    code1(files.con, files.json, files.raw, "", args["time-limit"], scoring_method, files.case_id, files.scenario_id, output_dir=args["case"])
    code1_time = time() - time_start

    if args["skip-solution2"]
        @warn("skip solution 2 solve")
        return
    end

    time_start = time()
    code2(files.con, files.json, files.raw, "", args["time-limit"], scoring_method, files.case_id, files.scenario_id, output_dir=args["case"])
    code2_time = time() - time_start

    evaluation_summary(args["case"], length(Distributed.workers()), code1_time, code2_time)

    if args["remove-solutions"]
        @warn("removing solution 2 files")
        remove_c2_solution_files(output_dir=args["case"])
    end
end


function parse_scopf_c2_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--case", "-c"
            help = "the case directory"
            required = true
        "--time-limit", "-t"
            help = "solution 1 runtime limit (seconds)"
            arg_type = Int
            default = 600000
        "--hard-time-limit"
            help = "solution 1 and 2 runtime limit (seconds)"
            arg_type = Int
            default = 600000
        "--skip-solution2"
            help = "skip computation of solution2 file"
            action = :store_true
        "--remove-solutions"
            help = "delete the solution files after competition"
            action = :store_true
    end

    return parse_args(s)
end


