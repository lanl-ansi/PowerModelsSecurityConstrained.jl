start_init = time()

include("distributed.jl")
add_procs()

println("cluster setup time: $(time() - start_init)")
println("")

start_pkg = time()

@everywhere include("goc_c2_huristic.jl")

println("package load time: $(time() - start_pkg)")
println("")

println("script startup time: $(time() - start_init)")
println("")

function MyJulia2(InFile1::String, InFile2::String, InFile3::String, InFile4::String, TimeLimitInSeconds::Int64, ScoringMethod::Int64, NetworkModel::String)
    println("running MyJulia2")
    println("  $(InFile1)")
    println("  $(InFile2)")
    println("  $(InFile3)")
    println("  $(InFile4)")
    println("  $(TimeLimitInSeconds)")
    println("  $(ScoringMethod)")
    println("  $(NetworkModel)")

    startup_time = time() - start_init
    remaining_time = trunc(Int, TimeLimitInSeconds-startup_time-60)
    println("")
    println("time limit: $(TimeLimitInSeconds)")
    println("startup time: $(startup_time)")
    println("remaining time: $(remaining_time)")

    code2(InFile1, InFile2, InFile3, InFile4, remaining_time, ScoringMethod, NetworkModel, "anonymous_scenario")
end
