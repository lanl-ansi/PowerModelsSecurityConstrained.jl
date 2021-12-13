#!/usr/bin/env julia

include("scopf-main.jl")
include("goc_c2_huristic.jl")

scopf_c2_main(parse_scopf_c2_commandline())
