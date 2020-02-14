#!/usr/bin/env julia

include("goc_challenge1_huristic.jl")
include("scopf-main.jl")

scopf_main(parse_scopf_commandline())
