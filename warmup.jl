#!/usr/bin/env julia

using Pkg
#Pkg.activate(".")
#Pkg.instantiate()
Pkg.status()

time_start = time()
println("start pre-compilation")

using Distributed
using SparseArrays

using JSON

using JuMP

using Ipopt

using InfrastructureModels
using PowerModels
using Memento

println("pre-compilation finished ($(time() - time_start))")
