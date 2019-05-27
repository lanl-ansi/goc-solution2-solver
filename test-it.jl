#!/usr/bin/env julia --project=.

using Distributed

InFile1="test-dataset/scenario_1/case.con"
InFile2="test-dataset/scenario_1/case.inl"
InFile3="test-dataset/scenario_1/case.raw"
InFile4="test-dataset/scenario_1/case.rop"
TimeLimitInSeconds=600
ScoringMethod=2
NetworkModel="IEEE 14"

include("solution2-solver.jl")

println("  $(InFile1)")
println("  $(InFile2)")
println("  $(InFile3)")
println("  $(InFile4)")
println("  $(TimeLimitInSeconds)")
println("  $(ScoringMethod)")
println("  $(NetworkModel)")

compute_solution2(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel, output_dir="test-dataset/scenario_1")
