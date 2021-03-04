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


ini_file = "../test/data/inputfiles.ini"
scenarios = ["scenario_01", "scenario_02"]
cases = [parse_goc_files(ini_file, scenario_id=sid) for sid in scenarios]

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
