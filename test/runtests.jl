using PowerModelsSecurityConstrained
using PowerModels

using Memento

using Test


Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsSecurityConstrained), "error")

ini_file = "../test/data/inputfiles.ini"
case_01 = parse_goc_files(ini_file, scenario_id="scenario_01")
case_02 = parse_goc_files(ini_file, scenario_id="scenario_02")

cases = [case_01, case_02]
networks = []

@testset "PowerModelsSecurityConstrained" begin

    include("contingency.jl")

end
