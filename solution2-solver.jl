@everywhere using Ipopt

include("code-2-lib/parsers.jl")
include("code-2-lib/lib.jl")

function compute_solution2(con_file::String, inl_file::String, raw_file::String, rop_file::String, time_limit::Int, scoring_method::Int, network_model::String; output_dir::String="", scenario_id::String="none")
    time_data_start = time()
    goc_data = parse_goc_files(con_file, inl_file, raw_file, rop_file, scenario_id=scenario_id)
    network = build_pm_model(goc_data)

    sol = read_solution1(network, output_dir=output_dir)
    PowerModels.update_data!(network, sol)

    check_network_solution(network)

    network_tmp = deepcopy(network)
    balance = compute_power_balance_deltas!(network_tmp)

    if balance.p_delta_abs_max > 0.01 || balance.q_delta_abs_max > 0.01
        error(LOGGER, "solution1 power balance requirements not satified (all power balance values should be below 0.01). $(balance)")
    end

    time_data = time() - time_data_start

    for (i,bus) in network["bus"]
        if haskey(bus, "evhi")
            bus["vmax"] = bus["evhi"]
        end
        if haskey(bus, "evlo")
            bus["vmin"] = bus["evlo"]
        end
    end

    for (i,branch) in network["branch"]
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        branch["g"] = g
        branch["b"] = b
        branch["tr"] = tr
        branch["ti"] = ti
    end


    ###### Prepare Solution 2 ######

    time_contingencies_start = time()

    gen_cont_total = length(network["gen_contingencies"])
    branch_cont_total = length(network["branch_contingencies"])

    processes = 1
    if Distributed.nprocs() > 1
        processes = Distributed.nprocs()-1 # save one for the master process
    end

    networks = [copy(network) for n in 1:processes]
    gens_per_proc = trunc(Int, ceil(gen_cont_total/processes))
    branch_per_proc = trunc(Int, ceil(branch_cont_total/processes))

    for p in 1:processes
        network_proc = networks[p]
        network_proc["gen_contingencies"] = network_proc["gen_contingencies"][(1+(p-1)*gens_per_proc):min(gen_cont_total,p*gens_per_proc)]
        network_proc["branch_contingencies"] = network_proc["branch_contingencies"][(1+(p-1)*branch_per_proc):min(branch_cont_total,p*branch_per_proc)]
    end

    for (i,network_proc) in enumerate(networks)
        info(LOGGER, "network $(i): $(length(network_proc["gen_contingencies"])), $(length(network_proc["branch_contingencies"]))")
    end

    result = pmap(solution2_solver, networks)

    contingency_solutions = Dict{String,Any}()
    for cont_sols in result
        for (id,sol) in cont_sols
            @assert !haskey(contingency_solutions, id)
            contingency_solutions[id] = sol
        end
    end

    correct_contingency_solutions!(network, contingency_solutions)
    write_solution2(network, contingency_solutions; output_dir=output_dir)

    time_contingencies = time() - time_contingencies_start
    info(LOGGER, "contingency eval time: $(time_contingencies)")
end


@everywhere function solution2_solver(network)
    bus_gens = gens_by_bus(network)

    network["delta"] = 0
    for (i,bus) in network["bus"]
        bus["vm_base"] = bus["vm"]
        bus["vm_start"] = bus["vm"]
        bus["va_start"] = bus["va"]
        bus["vm_fixed"] = length(bus_gens[i]) != 0

    end

    for (i,gen) in network["gen"]
        gen["pg_base"] = gen["pg"]
        gen["pg_start"] = gen["pg"]
        gen["qg_start"] = gen["qg"]
        gen["pg_fixed"] = false
        gen["qg_fixed"] = false
    end

    nlp_solver = JuMP.with_optimizer(Ipopt.Optimizer, tol=1e-6, print_level=0)

    contingency_solutions = Dict{String,Any}()

    for cont in network["gen_contingencies"]
        info(LOGGER, "working on: $(cont.label)")
        network_tmp = deepcopy(network)

        cont_gen = network_tmp["gen"]["$(cont.idx)"]
        cont_gen["contingency"] = true
        cont_gen["gen_status"] = 0
        pg_lost = cont_gen["pg"]

        result = run_fixpoint_pf_v2_2!(network_tmp, pg_lost, ACRPowerModel, nlp_solver, iteration_limit=5)

        result["solution"]["feasible"] = (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
        result["solution"]["cont_type"] = "gen"
        result["solution"]["cont_comp_id"] = cont.idx

        result["solution"]["gen"]["$(cont.idx)"]["pg"] = 0.0
        result["solution"]["gen"]["$(cont.idx)"]["qg"] = 0.0

        contingency_solutions[cont.label] = result["solution"]
        network_tmp["gen"]["$(cont.idx)"]["gen_status"] = 1
    end

    for cont in network["branch_contingencies"]
        info(LOGGER, "working on: $(cont.label)")
        network_tmp = deepcopy(network)
        network_tmp["branch"]["$(cont.idx)"]["br_status"] = 0

        result = run_fixpoint_pf_v2_2!(network_tmp, 0.0, ACRPowerModel, nlp_solver, iteration_limit=5)

        result["solution"]["feasible"] = (result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED)
        result["solution"]["cont_type"] = "branch"
        result["solution"]["cont_comp_id"] = cont.idx

        contingency_solutions[cont.label] = result["solution"]
        network_tmp["branch"]["$(cont.idx)"]["br_status"] = 1
    end

    return contingency_solutions
end
