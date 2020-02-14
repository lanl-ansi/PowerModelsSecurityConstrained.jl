start_init = time()

include("distributed.jl")
#add_procs() #can be restored after package registration

start_pkg = time()

#@everywhere using Pkg
#@everywhere Pkg.activate(".")

@everywhere using PowerModelsSecurityConstrained

@everywhere using JuMP
@everywhere using Ipopt
@everywhere using PowerModels

@everywhere using Memento
@everywhere const LOGGER = Memento.getlogger(PowerModelsSecurityConstrained)

include("goc_challenge1_huristic.jl")
include("second-stage-solution1.jl")

println("package load time: $(time() - start_pkg)")

println("script startup time: $(time() - start_init)")


function MyJulia2(InFile1::String, InFile2::String, InFile3::String, InFile4::String, TimeLimitInSeconds::Int64, ScoringMethod::Int64, NetworkModel::String)
    println("running MyJulia2")
    println("  $(InFile1)")
    println("  $(InFile2)")
    println("  $(InFile3)")
    println("  $(InFile4)")
    println("  $(TimeLimitInSeconds)")
    println("  $(ScoringMethod)")
    println("  $(NetworkModel)")

    compute_solution2_fast(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)

    compute_solution2(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)
end
