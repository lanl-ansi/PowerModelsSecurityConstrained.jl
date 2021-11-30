using PowerModelsSecurityConstrained
using PowerModels

using Memento

using Test

using Ipopt
using Cbc

Memento.setlevel!(Memento.getlogger(PowerModels.InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsSecurityConstrained), "error")


nlp_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)
#nlp_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6)
lp_solver = optimizer_with_attributes(Cbc.Optimizer, "logLevel"=>0)
#lp_solver = optimizer_with_attributes(Cbc.Optimizer)


c1_ini_file = "../test/data/c1/inputfiles.ini"
c1_scenarios = ["scenario_01", "scenario_02"]
c1_cases = [parse_c1_files(c1_ini_file, scenario_id=sid) for sid in c1_scenarios]

c1_networks = [build_c1_pm_model(case) for case in c1_cases]
c1_solutions = [read_c1_solution1(c1_networks[i], output_dir=dirname(case.files["raw"])) for (i,case) in enumerate(c1_cases)]

c1_case_infeasible = parse_c1_files(c1_ini_file, scenario_id="scenario_03")
c1_network_infeasible = build_c1_pm_model(c1_case_infeasible)


@testset "PowerModelsSecurityConstrained" begin

    include("common.jl")

    include("opf.jl")

    include("contingency.jl")

    include("contingency-stage.jl")

    include("scopf.jl")

end
