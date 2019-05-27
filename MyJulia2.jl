max_threads = 12

using Distributed
if distribute && Distributed.nprocs() <= 2
    threads = min(trunc(Int, Sys.CPU_THREADS*0.75), max_threads)
    println("threads: $(threads)/$(Sys.CPU_THREADS)")

    proc_ids = Distributed.addprocs(threads)
    println("process ids: $(proc_ids)")
end

#@everywhere using Pkg
#@everywhere Pkg.activate(".")

include("solution2-solver.jl")

function MyJulia2(InFile1::String, InFile2::String, InFile3::String, InFile4::String, TimeLimitInSeconds::Int64, ScoringMethod::Int64, NetworkModel::String)
    println("running MyJulia2")
    println("  $(InFile1)")
    println("  $(InFile2)")
    println("  $(InFile3)")
    println("  $(InFile4)")
    println("  $(TimeLimitInSeconds)")
    println("  $(ScoringMethod)")
    println("  $(NetworkModel)")

    compute_solution2(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)
end
