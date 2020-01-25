start_init = time()

include("distributed.jl")
add_procs()

start_pkg = time()

#@everywhere using Pkg
#@everywhere Pkg.activate(".")

@everywhere using PowerModelsSecurityConstrained

include("scopf-staged-solver-030.jl")

println("package load time: $(time() - start_pkg)")

println("script startup time: $(time() - start_init)")


function MyJulia1(InFile1::String, InFile2::String, InFile3::String, InFile4::String, TimeLimitInSeconds::Int64, ScoringMethod::Int64, NetworkModel::String)
    println("running MyJulia1")
    println("  $(InFile1)")
    println("  $(InFile2)")
    println("  $(InFile3)")
    println("  $(InFile4)")
    println("  $(TimeLimitInSeconds)")
    println("  $(ScoringMethod)")
    println("  $(NetworkModel)")

    startup_time = 60

    compute_solution1(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds-startup_time, ScoringMethod, NetworkModel)
end





