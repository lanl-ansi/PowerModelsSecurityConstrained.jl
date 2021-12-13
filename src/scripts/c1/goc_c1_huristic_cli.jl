#!/usr/bin/env julia

include("scopf-main.jl")
include("goc_c1_huristic.jl")

scopf_main(parse_scopf_commandline())
