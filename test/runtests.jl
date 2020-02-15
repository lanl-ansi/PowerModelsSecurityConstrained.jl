using PowerModelsSecurityConstrained
using PowerModels

using Memento

using Test

using Ipopt
using Cbc

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
solutions = [build_pm_solution(networks[i], joinpath(dirname(case.files["raw"]), "solution1.txt")) for (i,case) in enumerate(cases)]

@testset "PowerModelsSecurityConstrained" begin

    include("opf.jl")

    include("contingency.jl")

    include("security-stage.jl")

    include("scopf.jl")

end
