#!/usr/bin/env julia --project=.

InFile1="../../test/data/scenario_01/case.con"
InFile2="../../test/data/scenario_01/case.inl"
InFile3="../../test/data/scenario_01/case.raw"
InFile4="../../test/data/scenario_01/case.rop"
TimeLimitInSeconds=600
ScoringMethod=2
NetworkModel="IEEE 14"

include("MyJulia1.jl")
MyJulia1(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)

include("MyJulia2.jl")
MyJulia2(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)
