using PowerModelsSecurityConstrained
using PowerModels

using Memento

using Test

using Ipopt
using Cbc

Memento.setlevel!(Memento.getlogger(PowerModels.InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsSecurityConstrained), "error")


nlp_solver = with_optimizer(Ipopt.Optimizer, tol=1e-6, print_level=0)
#nlp_solver = with_optimizer(Ipopt.Optimizer, tol=1e-6)
lp_solver = with_optimizer(Cbc.Optimizer, logLevel=0)
#lp_solver = with_optimizer(Cbc.Optimizer)


ini_file = "../test/data/inputfiles.ini"
case_01 = parse_goc_files(ini_file, scenario_id="scenario_01")
case_02 = parse_goc_files(ini_file, scenario_id="scenario_02")
cases = [case_01, case_02]

networks = [build_pm_model(case) for case in cases]
solutions = [read_solution1(networks[i], output_dir=dirname(case.files["raw"])) for (i,case) in enumerate(cases)]

case_infeasible = parse_goc_files(ini_file, scenario_id="scenario_03")
network_infeasible = build_pm_model(case_infeasible)

@testset "PowerModelsSecurityConstrained" begin

    include("common.jl")

    include("opf.jl")

    include("contingency.jl")

    include("contingency-stage.jl")

    include("scopf.jl")

end
