@everywhere using JuMP

@everywhere using PowerModels

@everywhere using InfrastructureModels

@everywhere using Memento

@everywhere const LOGGER = getlogger("GOC")

@everywhere const vm_eq_tol = 1e-4

@everywhere import Statistics: mean


@everywhere function build_pm_model(goc_data)
    scenario = goc_data.scenario
    network = goc_data.network

    ##### General Helpers #####
    gen_lookup = Dict(tuple(gen["source_id"][2], strip(gen["source_id"][3])) => gen for (i,gen) in network["gen"])

    branch_lookup = Dict()
    for (i,branch) in network["branch"]
        if !branch["transformer"]
            branch_id = tuple(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][4]))
        else
            branch_id = tuple(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][5]))
            @assert branch["source_id"][4] == 0
            @assert branch["source_id"][6] == 0
        end
        branch_lookup[branch_id] = branch
    end



    ##### Link Generator Cost Data #####

    @assert network["per_unit"]
    mva_base = network["baseMVA"]

    dispatch_tbl_lookup = Dict()
    for dispatch_tbl in goc_data.cost["disptbl"]
        dispatch_tbl_lookup[dispatch_tbl["ctbl"]] = dispatch_tbl
    end

    cost_tbl_lookup = Dict()
    for cost_tbl in goc_data.cost["ctbl"]
        cost_tbl_lookup[cost_tbl["ltbl"]] = cost_tbl
    end

    gen_cost_models = Dict()
    for gen_dispatch in goc_data.cost["gen"]
        gen_id = (gen_dispatch["bus"], strip(gen_dispatch["genid"]))
        dispatch_tbl = dispatch_tbl_lookup[gen_dispatch["disptbl"]]
        cost_tbl = cost_tbl_lookup[dispatch_tbl["ctbl"]]

        gen_cost_models[gen_id] = cost_tbl
    end

    if length(gen_cost_models) != length(network["gen"])
        error(LOGGER, "cost model data missing, network has $(length(network["gen"])) generators, the cost model has $(length(gen_cost_models)) generators")
    end

    for (gen_id, cost_model) in gen_cost_models
        pm_gen = gen_lookup[gen_id]
        pm_gen["model"] = 1
        pm_gen["model_lable"] = cost_model["label"]
        pm_gen["ncost"] = length(cost_model["points"])

        #println(cost_model["points"])
        point_list = Float64[]
        for point in cost_model["points"]
            push!(point_list, point.x/mva_base)
            push!(point_list, point.y)
        end
        pm_gen["cost"] = point_list
    end



    ##### Link Generator Participation Data #####

    if length(goc_data.response) != length(network["gen"])
        error(LOGGER, "generator response model data missing, network has $(length(network["gen"])) generators, the response model has $(length(goc_data.response)) generators")
    end

    for gen_response in goc_data.response
        gen_id = (gen_response["i"], strip(gen_response["id"]))

        pm_gen = gen_lookup[gen_id]

        pm_gen["alpha"] = gen_response["r"]
    end



    ##### Flexible Shunt Data #####

    for (i,shunt) in network["shunt"]
        if shunt["source_id"][1] == "switched shunt"
            @assert shunt["source_id"][3] == 0
            @assert shunt["gs"] == 0.0
            shunt["dispatchable"] = true

            bmin = 0.0
            bmax = 0.0
            for (n_name,b_name) in [("n1","b1"),("n2","b2"),("n3","b3"),("n4","b4"),("n5","b5"),("n6","b6"),("n7","b7"),("n8","b8")]
                if shunt[b_name] <= 0.0
                    bmin += shunt[n_name]*shunt[b_name]
                else
                    bmax += shunt[n_name]*shunt[b_name]
                end
            end
            shunt["bmin"] = bmin/mva_base
            shunt["bmax"] = bmax/mva_base
        else
            shunt["dispatchable"] = false
        end
    end



    ##### Add Contingency Lists #####

    generator_ids = []
    branch_ids = []

    for (i,cont) in enumerate(goc_data.contingencies)
        if cont["component"] == "branch"
            branch_id = (cont["i"], cont["j"], cont["ckt"])
            pm_branch = branch_lookup[branch_id]
            push!(branch_ids, (idx=pm_branch["index"], label=cont["label"], type="branch"))

        elseif cont["component"] == "generator"
            gen_id = (cont["i"], cont["id"])
            pm_gen = gen_lookup[gen_id]
            push!(generator_ids, (idx=pm_gen["index"], label=cont["label"], type="gen"))

        else
            error(LOGGER, "unrecognized contingency component type $(cont["component"]) at contingency $(i)")
        end
    end

    network["branch_contingencies"] = branch_ids
    network["gen_contingencies"] = generator_ids

    network["branch_contingencies_active"] = []
    network["gen_contingencies_active"] = []



    ##### Fix Broken Data #####

    PowerModels.correct_cost_functions!(network)

    # FYI, this breaks output API
    #PowerModels.propagate_topology_status!(network)

    for (i,shunt) in network["shunt"]
        # test checks if a "switched shunt" in the orginal data model
        if shunt["dispatchable"]
            if shunt["bs"] < shunt["bmin"]
                warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmin"])")
                shunt["bs"] = shunt["bmin"]
            end
            if shunt["bs"] > shunt["bmax"]
                warn(LOGGER, "update bs on shunt $(i) to be in bounds $(shunt["bs"]) -> $(shunt["bmax"])")
                shunt["bs"] = shunt["bmax"]
            end
        end
    end

    return network
end


function build_pm_opf_model(goc_data)
    scenario = goc_data.scenario
    network = goc_data.network

    ##### General Helpers #####

    gen_lookup = Dict(tuple(gen["source_id"][2], strip(gen["source_id"][3])) => gen for (i,gen) in network["gen"])

    branch_lookup = Dict()
    for (i,branch) in network["branch"]
        if !branch["transformer"]
            branch_id = tuple(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][4]))
        else
            branch_id = tuple(branch["source_id"][2], branch["source_id"][3], strip(branch["source_id"][5]))
            @assert branch["source_id"][4] == 0
            @assert branch["source_id"][6] == 0
        end
        branch_lookup[branch_id] = branch
    end



    ##### Link Generator Cost Data #####

    @assert network["per_unit"]
    mva_base = network["baseMVA"]

    dispatch_tbl_lookup = Dict()
    for dispatch_tbl in goc_data.cost["disptbl"]
        dispatch_tbl_lookup[dispatch_tbl["ctbl"]] = dispatch_tbl
    end

    cost_tbl_lookup = Dict()
    for cost_tbl in goc_data.cost["ctbl"]
        cost_tbl_lookup[cost_tbl["ltbl"]] = cost_tbl
    end

    gen_cost_models = Dict()
    for gen_dispatch in goc_data.cost["gen"]
        gen_id = (gen_dispatch["bus"], strip(gen_dispatch["genid"]))
        dispatch_tbl = dispatch_tbl_lookup[gen_dispatch["disptbl"]]
        cost_tbl = cost_tbl_lookup[dispatch_tbl["ctbl"]]

        gen_cost_models[gen_id] = cost_tbl
    end

    if length(gen_cost_models) != length(network["gen"])
        error(LOGGER, "cost model data missing, network has $(length(network["gen"])) generators, the cost model has $(length(gen_cost_models)) generators")
    end

    for (gen_id, cost_model) in gen_cost_models
        pm_gen = gen_lookup[gen_id]
        pm_gen["model"] = 1
        pm_gen["model_lable"] = cost_model["label"]
        pm_gen["ncost"] = length(cost_model["points"])

        #println(cost_model["points"])
        point_list = Float64[]
        for point in cost_model["points"]
            push!(point_list, point.x/mva_base)
            push!(point_list, point.y)
        end
        pm_gen["cost"] = point_list
    end

    PowerModels.correct_cost_functions!(network)

    return network
end



@everywhere function read_solution1(network; output_dir="", state_file="solution1.txt")
    if length(output_dir) > 0
        solution1_path = joinpath(output_dir, state_file)
    else
        solution1_path = state_file
    end

    return build_pm_solution(network, solution1_path)
end

@everywhere function build_pm_solution(network, goc_sol_file::String)
    info(LOGGER, "loading solution file: $(goc_sol_file)")
    goc_sol = parse_solution1_file(goc_sol_file)

    info(LOGGER, "converting GOC solution to PowerModels solution")
    pm_sol = build_pm_solution(network, goc_sol)

    return pm_sol
end

@everywhere function build_pm_solution(network, goc_sol)
    bus_lookup = Dict(parse(Int, bus["source_id"][2]) => bus for (i,bus) in network["bus"])
    gen_lookup = Dict((gen["source_id"][2], strip(gen["source_id"][3])) => gen for (i,gen) in network["gen"])
    shunt_lookup = Dict{Int,Any}()
    for (i,shunt) in network["shunt"]
        if shunt["source_id"][1] == "switched shunt"
            @assert shunt["source_id"][3] == 0
            shunt_lookup[shunt["source_id"][2]] = shunt
        end
    end

    base_mva = network["baseMVA"]

    bus_data = Dict{String,Any}()
    shunt_data = Dict{String,Any}()
    for bus_sol in goc_sol.bus
        pm_bus = bus_lookup[bus_sol.bus]
        bus_data["$(pm_bus["index"])"] = Dict(
            "vm" => bus_sol.vm,
            "va" => deg2rad(bus_sol.va)
        )

        if haskey(shunt_lookup, bus_sol.bus)
            pm_shunt = shunt_lookup[bus_sol.bus]
            shunt_data["$(pm_shunt["index"])"] = Dict(
                "gs" => 0.0,
                "bs" => bus_sol.bcs/base_mva
            )
        else
            @assert bus_sol.bcs == 0.0
        end
    end

    gen_data = Dict{String,Any}()
    for gen_sol in goc_sol.gen
        pm_gen = gen_lookup[(gen_sol.bus, gen_sol.id)]
        gen_data["$(pm_gen["index"])"] = Dict(
            "pg" => gen_sol.pg/base_mva,
            "qg" => gen_sol.qg/base_mva
        )
    end

    solution = Dict(
        "per_unit" => true,
        "bus" => bus_data,
        "shunt" => shunt_data,
        "gen" => gen_data
    )

    return solution
end

@everywhere function gens_by_bus(network)
    bus_gens = Dict(i => Any[] for (i,bus) in network["bus"])
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            push!(bus_gens["$(gen["gen_bus"])"], gen)
        end
    end
    return bus_gens
end



@everywhere function run_fixpoint_pf_v2_2!(network, pg_lost, model_constructor, solver; iteration_limit=typemax(Int64))
    time_start = time()

    delta = apply_pg_response!(network, pg_lost)
    #delta = 0.0

    debug(LOGGER, "pg lost: $(pg_lost)")
    debug(LOGGER, "delta: $(network["delta"])")
    debug(LOGGER, "pre-solve time: $(time() - time_start)")

    base_solution = extract_solution(network)
    base_solution["delta"] = delta
    final_result = Dict(
        "termination_status" => LOCALLY_SOLVED,
        "solution" => base_solution
    )

    bus_gens = gens_by_bus(network)

    #network["delta_start"] = network["delta"]
    for (i,gen) in network["gen"]
        gen["qg_fixed"] = false
        gen["pg_start"] = gen["pg"]
        if isapprox(gen["qmin"],gen["qmax"])
            gen["qg_fixed"] = true
            gen["qg"] = gen["qmin"]
        end
        gen["qg_start"] = gen["qg"]
    end

    for (i,bus) in network["bus"]
        active_gens = [gen for gen in bus_gens[i] if !gen["qg_fixed"]]
        if length(active_gens) == 0
            bus["vm_fixed"] = false
        else
            bus["vm_fixed"] = true
        end
        #bus["vm_start"] = bus["vm"]
        #bus["va_start"] = bus["va"]
        bus["vr_start"] = bus["vm"]*cos(bus["va"])
        bus["vi_start"] = bus["vm"]*sin(bus["va"])
    end

    time_start = time()
    result = run_fixed_pf_nbf_rect2(network, model_constructor, solver)
    debug(LOGGER, "pf solve time: $(time() - time_start)")
    if result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
        correct_qg!(network, result["solution"], bus_gens=bus_gens)
        PowerModels.update_data!(network, result["solution"])
        final_result = result
    else
        warn(LOGGER, "contingency pf solver FAILED with status $(result["termination_status"]) on iteration 0")
        return final_result
    end






    pg_switched = true
    qg_switched = true
    vm_switched = true

    iteration = 1
    deltas = [result["solution"]["delta"]]
    while (pg_switched || qg_switched || vm_switched) && iteration <= iteration_limit
        debug(LOGGER, "obj: $(result["objective"])")
        debug(LOGGER, "delta: $(result["solution"]["delta"])")
        pg_switched = false
        qg_switched = false
        vm_switched = false

        for (i,gen) in network["gen"]
            pg = gen["pg_base"] + network["delta"]*gen["alpha"]

            if gen["pg_fixed"]
                if !isapprox(gen["pmax"], gen["pmin"]) && pg < gen["pmax"] && pg > gen["pmin"]
                    gen["pg"] = pg
                    gen["pg_fixed"] = false
                    pg_switched = true
                    #info(LOGGER, "unfix pg on gen $(i)")
                end
            else
                if pg >= gen["pmax"]
                    gen["pg"] = gen["pmax"]
                    gen["pg_fixed"] = true
                    pg_switched = true
                    #info(LOGGER, "fix pg to ub on gen $(i)")
                elseif gen["pg"] <= gen["pmin"]
                    gen["pg"] = gen["pmin"]
                    gen["pg_fixed"] = true
                    pg_switched = true
                    #info(LOGGER, "fix pg to lb on gen $(i)")
                end
            end
        end

        for (i,bus) in network["bus"]
            if length(bus_gens[i]) > 0
                qg = sum(gen["qg"] for gen in bus_gens[i])
                qmin = sum(gen["qmin"] for gen in bus_gens[i])
                qmax = sum(gen["qmax"] for gen in bus_gens[i])

                if isapprox(qmin,qmax)
                    @assert !bus["vm_fixed"]
                    for gen in bus_gens[i]
                        @assert gen["qg_fixed"]
                        @assert isapprox(gen["qg"],gen["qmin"])
                    end
                elseif bus["vm_fixed"]
                    if qg >= qmax
                        bus["vm_fixed"] = false
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg"] = gen["qmax"]
                            gen["qg_fixed"] = true
                        end
                        #info(LOGGER, "fix qg to ub on bus $(i)")
                    end

                    if qg <= qmin
                        bus["vm_fixed"] = false
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg"] = gen["qmin"]
                            gen["qg_fixed"] = true
                        end
                        #info(LOGGER, "fix qg to lb on bus $(i)")
                    end
                else
                    if qg < qmax && qg > qmin
                        bus["vm_fixed"] = true

                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg_fixed"] = false
                            gen["qg_start"] = gen["qg"]
                        end
                        #info(LOGGER, "fix vm to on bus $(i)")
                    end
                    if qg >= qmax && bus["vm"] > bus["vm_base"]
                        bus["vm_fixed"] = true
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg_fixed"] = false
                        end
                        #info(LOGGER, "fix vm to on bus $(i)")
                    end
                    if qg <= qmin && bus["vm"] < bus["vm_base"]
                        bus["vm_fixed"] = true
                        vm_switched = true
                        for gen in bus_gens[i]
                            gen["qg_fixed"] = false
                        end
                        #info(LOGGER, "fix vm to on bus $(i)")
                    end
                end
            end
        end


        for (i,gen) in network["gen"]
            gen["pg_start"] = gen["pg"]
            gen["qg_start"] = gen["qg"]
        end

        for (i,bus) in network["bus"]
            bus["vm_start"] = bus["vm"]
            bus["va_start"] = bus["va"]
        end


        if pg_switched || qg_switched || vm_switched
            debug(LOGGER, "bus or gen swtiched: $iteration")
            time_start = time()
            result = run_fixed_pf_nbf_rect2(network, model_constructor, solver)
            debug(LOGGER, "pf solve time: $(time() - time_start)")
            if result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
                correct_qg!(network, result["solution"], bus_gens=bus_gens)
                PowerModels.update_data!(network, result["solution"])
                final_result = result
            else
                warn(LOGGER, "contingency pf solver FAILED with status $(result["termination_status"]) on iteration 0")
                return final_result
            end

            push!(deltas, result["solution"]["delta"])
            iteration += 1
            if iteration >= iteration_limit
                warn(LOGGER, "hit iteration limit")
            end
            if length(deltas) > 3 && isapprox(deltas[end-2], deltas[end])
                warn(LOGGER, "cycle detected, stopping")
                break
            end
        end
    end

    return final_result
end



@everywhere function apply_pg_response!(network, pg_delta)
    for (i,gen) in network["gen"]
        gen["pg_fixed"] = false
    end

    pg_total = 0.0
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            pg_total += gen["pg"]
        end
    end

    pg_target = pg_total + pg_delta
    #info(LOGGER, "total gen:  $(pg_total)")
    #info(LOGGER, "target gen: $(pg_target)")


    delta_est = 0.0
    while !isapprox(pg_total, pg_target)
        alpha_total = 0.0
        for (i,gen) in network["gen"]
            if gen["gen_status"] != 0 && !gen["pg_fixed"]
                alpha_total += gen["alpha"]
            end
        end
        #info(LOGGER, "alpha total: $(alpha_total)")

        if isapprox(alpha_total, 0.0) && !isapprox(pg_total, pg_target)
            warn(LOGGER, "insufficient generator response to meet demand, remaining pg $(pg_total - pg_target), remaining alpha $(alpha_total)")
            break
        end

        delta_est += pg_delta/alpha_total
        #info(LOGGER, "detla: $(delta_est)")

        for (i,gen) in network["gen"]
            if gen["gen_status"] != 0
                pg_cont = gen["pg_base"] + delta_est*gen["alpha"]

                if pg_cont <= gen["pmin"]
                    gen["pg"] = gen["pmin"]
                    if !gen["pg_fixed"]
                        gen["pg_fixed"] = true
                        #info(LOGGER, "gen $(i) hit lb $(gen["pmin"]) with target value of $(pg_cont)")
                    end
                elseif pg_cont >= gen["pmax"]
                    gen["pg"] = gen["pmax"]
                    if !gen["pg_fixed"]
                        gen["pg_fixed"] = true
                        #info(LOGGER, "gen $(i) hit ub $(gen["pmax"]) with target value of $(pg_cont)")
                    end
                else
                    gen["pg"] = pg_cont
                end
            end
        end

        pg_total = 0.0
        for (i,gen) in network["gen"]
            if gen["gen_status"] != 0
                pg_total += gen["pg"]
            end
        end
        #pg_comp = comp_pg_response_total(network, delta_est)
        #info(LOGGER, "detla: $(delta_est)")
        #info(LOGGER, "total gen comp $(pg_comp) - gen inc $(pg_total)")
        #info(LOGGER, "total gen $(pg_total) - target gen $(pg_target)")

        pg_delta = pg_target - pg_total
    end

    network["delta"] = delta_est
    return delta_est
end



@everywhere function extract_solution(network; branch_flow=false)
    sol = Dict{String,Any}()

    sol["bus"] = Dict{String,Any}()
    for (i,bus) in network["bus"]
        bus_dict = Dict{String,Any}()
        bus_dict["va"] = get(bus, "va", 0.0)
        bus_dict["vm"] = get(bus, "vm", 1.0)
        sol["bus"][i] = bus_dict
    end

    sol["shunt"] = Dict{String,Any}()
    for (i,shunt) in network["shunt"]
        shunt_dict = Dict{String,Any}()
        shunt_dict["gs"] = get(shunt, "gs", 0.0)
        shunt_dict["bs"] = get(shunt, "bs", 0.0)
        sol["shunt"][i] = shunt_dict
    end

    sol["gen"] = Dict{String,Any}()
    for (i,gen) in network["gen"]
        gen_dict = Dict{String,Any}()
        gen_dict["pg"] = get(gen, "pg", 0.0)
        gen_dict["qg"] = get(gen, "qg", 0.0)
        sol["gen"][i] = gen_dict
    end

    if branch_flow
        sol["branch"] = Dict{String,Any}()
        for (i,branch) in network["branch"]
            branch_dict = Dict{String,Any}()
            branch_dict["pf"] = get(branch, "pf", 0.0)
            branch_dict["qf"] = get(branch, "qf", 0.0)
            branch_dict["pt"] = get(branch, "pt", 0.0)
            branch_dict["qt"] = get(branch, "qt", 0.0)
            sol["branch"][i] = branch_dict
        end
    end

    return sol
end


""

@everywhere function run_fixed_pf_nbf_rect2(file, model_constructor, solver; kwargs...)
    return run_model(file, model_constructor, solver, post_fixed_pf_nbf_rect2; solution_builder = solution_second_stage!, kwargs...)
end


""

@everywhere function post_fixed_pf_nbf_rect2(pm::GenericPowerModel)
    start_time = time()
    PowerModels.variable_voltage(pm, bounded=false)
    PowerModels.variable_active_generation(pm, bounded=false)
    PowerModels.variable_reactive_generation(pm, bounded=false)
    #PowerModels.variable_branch_flow(pm, bounded=false)
    #PowerModels.variable_dcline_flow(pm, bounded=false)

    #PowerModels.variable_branch_flow(pm, bounded=false)

    # TODO set bounds bounds on alpha and total gen capacity
    var(pm)[:delta] = @variable(pm.model, delta, base_name="delta", start=0.0)
    #var(pm)[:delta] = @variable(pm.model, delta, base_name="delta", start=ref(pm, :delta_start))
    #Memento.info(LOGGER, "post variable time: $(time() - start_time)")

    start_time = time()

    vr = var(pm, :vr)
    vi = var(pm, :vi)

    PowerModels.constraint_model_voltage(pm)

    for (i,bus) in ref(pm, :bus)
        if bus["vm_fixed"]
            @constraint(pm.model, vr[i]^2 + vi[i]^2 == bus["vm_base"]^2)
        end
    end

    for i in ids(pm, :ref_buses)
        PowerModels.constraint_theta_ref(pm, i)
    end
    #Memento.info(LOGGER, "misc constraints time: $(time() - start_time)")


    start_time = time()
    p = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    q = Dict{Tuple{Int64,Int64,Int64},GenericQuadExpr{Float64,VariableRef}}()
    for (i,branch) in ref(pm, :branch)
        #PowerModels.constraint_ohms_yt_from(pm, i)
        #PowerModels.constraint_ohms_yt_to(pm, i)

        f_bus_id = branch["f_bus"]
        t_bus_id = branch["t_bus"]
        f_idx = (i, f_bus_id, t_bus_id)
        t_idx = (i, t_bus_id, f_bus_id)

        f_bus = ref(pm, :bus, f_bus_id)
        t_bus = ref(pm, :bus, t_bus_id)

        #g, b = PowerModels.calc_branch_y(branch)
        #tr, ti = PowerModels.calc_branch_t(branch)
        g = branch["g"]
        b = branch["b"]
        tr = branch["tr"]
        ti = branch["ti"]

        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        #p_fr  = var(pm, :p, f_idx)
        #q_fr  = var(pm, :q, f_idx)
        #p_to  = var(pm, :p, t_idx)
        #q_to  = var(pm, :q, t_idx)
        vr_fr = vr[f_bus_id] #var(pm, :vr, f_bus_id)
        vr_to = vr[t_bus_id] #var(pm, :vr, t_bus_id)
        vi_fr = vi[f_bus_id] #var(pm, :vi, f_bus_id)
        vi_to = vi[t_bus_id] #var(pm, :vi, t_bus_id)

        p[f_idx] = (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        p[t_idx] = (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
    end
    #Memento.info(LOGGER, "flow expr time: $(time() - start_time)")


    start_time = time()
    pg = Dict{Int,Any}()
    qg = Dict{Int,Any}()
    for (i,gen) in ref(pm, :gen)
        if gen["pg_fixed"]
            pg[i] = gen["pg"]
            @constraint(pm.model, var(pm, :pg, i) == gen["pg"])
        else
            pg[i] = gen["pg_base"] + gen["alpha"]*delta
            @constraint(pm.model, var(pm, :pg, i) == gen["pg_base"] + gen["alpha"]*delta)
        end

        if gen["qg_fixed"]
            qg[i] = gen["qg"]
            @constraint(pm.model, var(pm, :qg, i) == gen["qg"])
        else
            qg[i] = var(pm, :qg, i)
        end
    end
    #Memento.info(LOGGER, "gen expr time: $(time() - start_time)")


    start_time = time()
    for (i,bus) in ref(pm, :bus)
        #PowerModels.constraint_kcl_shunt(pm, i)

        bus_arcs = ref(pm, :bus_arcs, i)
        bus_arcs_dc = ref(pm, :bus_arcs_dc, i)
        bus_gens = ref(pm, :bus_gens, i)
        bus_loads = ref(pm, :bus_loads, i)
        bus_shunts = ref(pm, :bus_shunts, i)

        bus_pd = Dict(k => ref(pm, :load, k, "pd") for k in bus_loads)
        bus_qd = Dict(k => ref(pm, :load, k, "qd") for k in bus_loads)

        bus_gs = Dict(k => ref(pm, :shunt, k, "gs") for k in bus_shunts)
        bus_bs = Dict(k => ref(pm, :shunt, k, "bs") for k in bus_shunts)

        #p = var(pm, :p)
        #q = var(pm, :q)
        #pg = var(pm, :pg)
        #qg = var(pm, :qg)

        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*(vr[i]^2 + vi[i]^2))
        @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*(vr[i]^2 + vi[i]^2))
    end
    #Memento.info(LOGGER, "power balance constraint time: $(time() - start_time)")

end


""

@everywhere function solution_second_stage!(pm::GenericPowerModel, sol::Dict{String,Any})
    #start_time = time()
    PowerModels.add_setpoint_bus_voltage!(sol, pm)
    #Memento.info(LOGGER, "voltage solution time: $(time() - start_time)")

    #start_time = time()
    PowerModels.add_setpoint_generator_power!(sol, pm)
    #Memento.info(LOGGER, "generator solution time: $(time() - start_time)")

    #start_time = time()
    PowerModels.add_setpoint_branch_flow!(sol, pm)
    #Memento.info(LOGGER, "branch solution time: $(time() - start_time)")

    sol["delta"] = JuMP.value(var(pm, :delta))

    #start_time = time()
    PowerModels.add_setpoint_fixed!(sol, pm, "shunt", "bs", default_value = (item) -> item["bs"])
    #Memento.info(LOGGER, "shunt solution time: $(time() - start_time)")

    #PowerModels.print_summary(sol)
end


"fixes solution degeneracy issues when qg is a free variable, as is the case in PowerFlow"

@everywhere function correct_qg!(network, solution; bus_gens=gens_by_bus(network))
    for (i,gens) in bus_gens
        if length(gens) > 1
            gen_ids = [gen["index"] for gen in gens]
            qgs = [solution["gen"]["$(j)"]["qg"] for j in gen_ids]
            if !isapprox(abs(sum(qgs)), sum(abs.(qgs)))
                #info(LOGGER, "$(i) - $(gen_ids) - $(qgs) - output requires correction!")
                qg_total = sum(qgs)

                qg_remaining = sum(qgs)
                qg_assignment = Dict(j => 0.0 for j in gen_ids)
                for (i,gen) in enumerate(gens)
                    gen_qg = qg_remaining
                    gen_qg = max(gen_qg, gen["qmin"])
                    gen_qg = min(gen_qg, gen["qmax"])
                    qg_assignment[gen["index"]] = gen_qg
                    qg_remaining = qg_remaining - gen_qg
                    if i == length(gens) && abs(qg_remaining) > 0.0
                        qg_assignment[gen["index"]] = gen_qg + qg_remaining
                    end
                end
                #info(LOGGER, "$(qg_assignment)")
                for (j,qg) in qg_assignment
                    solution["gen"]["$(j)"]["qg"] = qg
                end

                sol_qg_total = sum(solution["gen"]["$(j)"]["qg"] for j in gen_ids)
                #info(LOGGER, "$(qg_total) - $(sol_qg_total)")
                @assert isapprox(qg_total, sol_qg_total)

                #info(LOGGER, "updated to $([solution["gen"]["$(j)"]["qg"] for j in gen_ids])")
            end
        end
    end
end



"build a static ordering of all contigencies"

@everywhere function contingency_order(pm_network)
    gen_cont_order = sort(pm_network["gen_contingencies"], by=(x) -> x.label)
    branch_cont_order = sort(pm_network["branch_contingencies"], by=(x) -> x.label)

    gen_cont_total = length(gen_cont_order)
    branch_cont_total = length(branch_cont_order)

    gen_rate = 1.0
    branch_rate = 1.0
    steps = 1

    if gen_cont_total == 0 && branch_cont_total == 0
        # defaults are good
    elseif gen_cont_total == 0 && branch_cont_total != 0
        steps = branch_cont_total
    elseif gen_cont_total != 0 < branch_cont_total == 0
        steps = gen_cont_total
    elseif gen_cont_total == branch_cont_total
        steps = branch_cont_total
    elseif gen_cont_total < branch_cont_total
        gen_rate = 1.0
        branch_rate = branch_cont_total/gen_cont_total
        steps = gen_cont_total
    elseif gen_cont_total > branch_cont_total
        gen_rate = gen_cont_total/branch_cont_total
        branch_rate = 1.0 
        steps = branch_cont_total
    end

    cont_order = []
    gen_cont_start = 1
    branch_cont_start = 1
    for s in 1:steps
        gen_cont_end = min(gen_cont_total, trunc(Int,ceil(s*gen_rate)))
        #println(gen_cont_start:gen_cont_end)
        for j in gen_cont_start:gen_cont_end
            push!(cont_order, gen_cont_order[j])
        end
        gen_cont_start = gen_cont_end+1

        branch_cont_end = min(branch_cont_total, trunc(Int,ceil(s*branch_rate)))
        #println("$(s) - $(branch_cont_start:branch_cont_end)")
        for j in branch_cont_start:branch_cont_end
            push!(cont_order, branch_cont_order[j])
        end
        branch_cont_start = branch_cont_end+1
    end

    @assert(length(cont_order) == gen_cont_total + branch_cont_total)

    return cont_order
end



"checks feasibility criteria of contingencies, corrects when possible"

@everywhere function correct_contingency_solution!(network, cont_sol; bus_gens = gens_by_bus(network))
    label = cont_sol["label"]
    vm_changes = [0.0]
    for (i,bus) in cont_sol["bus"]
        nw_bus = network["bus"][i]

        if nw_bus["bus_type"] != 4
            if length(bus_gens[i]) > 0
                qg = sum(cont_sol["gen"]["$(gen["index"])"]["qg"] for gen in bus_gens[i])
                qmin = sum(gen["qmin"] for gen in bus_gens[i])
                qmax = sum(gen["qmax"] for gen in bus_gens[i])

                if !isapprox(abs(qmin - qmax), 0.0)
                    if qg >= qmax && bus["vm"] > nw_bus["vm"]
                        #println(qmin, " ", qg, " ", qmax)
                        #warn(LOGGER, "update qg $(qmin) $(qg) $(qmax)")
                        warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match set-point $(bus["vm"]) -> $(nw_bus["vm"]) due to qg upper bound and vm direction")
                        push!(vm_changes, abs(bus["vm"] - nw_bus["vm"]))
                        bus["vm"] = nw_bus["vm"]
                    end

                    if qg <= qmin && bus["vm"] < nw_bus["vm"]
                        #println(qmin, " ", qg, " ", qmax)
                        #warn(LOGGER, "update qg $(qmin) $(qg) $(qmax)")
                        warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match set-point $(bus["vm"]) -> $(nw_bus["vm"]) due to qg lower bound and vm direction")
                        push!(vm_changes, abs(bus["vm"] - nw_bus["vm"]))
                        bus["vm"] = nw_bus["vm"]
                    end
                end
            end

            if bus["vm"] > nw_bus["vmax"]
                warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match ub $(bus["vm"]) -> $(nw_bus["vmax"]) due to out of bounds")
                push!(vm_changes, bus["vm"] - nw_bus["vmax"])
                bus["vm"] = nw_bus["vmax"]
            end

            if bus["vm"] < nw_bus["vmin"]
                warn(LOGGER, "update vm on bus $(i) in contingency $(label) to match lb $(bus["vm"]) -> $(nw_bus["vmin"]) due to out of bounds")
                push!(vm_changes, nw_bus["vmin"] - bus["vm"])
                bus["vm"] = nw_bus["vmin"]
            end
        else
            bus["vm"] = 0.0
            bus["va"] = 0.0
        end
    end


    bs_changes = [0.0]
    for (i,shunt) in cont_sol["shunt"]
        nw_shunt = network["shunt"][i]
        if haskey(nw_shunt, "dispatchable") && nw_shunt["dispatchable"]
            @assert nw_shunt["gs"] == 0.0
            @assert haskey(nw_shunt, "bmin") && haskey(nw_shunt, "bmax")
            if shunt["bs"] > nw_shunt["bmax"]
                warn(LOGGER, "update bs on shunt $(i) in contingency $(label) to be in bounds $(shunt["bs"]) -> $(nw_shunt["bmax"])")
                push!(bs_changes, shunt["bs"] - nw_shunt["bmax"])
                shunt["bs"] = nw_shunt["bmax"]
            end
            if shunt["bs"] < nw_shunt["bmin"]
                warn(LOGGER, "update bs on shunt $(i) in contingency $(label) to be in bounds $(shunt["bs"]) -> $(nw_shunt["bmin"])")
                push!(bs_changes, nw_shunt["bmin"] - shunt["bs"])
                shunt["bs"] = nw_shunt["bmin"]
            end
        end
    end

    gen_id = -1
    if cont_sol["cont_type"] == "gen"
        gen_id = cont_sol["cont_comp_id"]
    end

    pg_changes = [0.0]
    qg_changes = [0.0]
    delta = cont_sol["delta"]
    for (i,gen) in cont_sol["gen"]
        nw_gen = network["gen"][i]

        if !(nw_gen["gen_status"] == 0 || (gen_id >= 0 && nw_gen["index"] == gen_id))
            bus_id = nw_gen["gen_bus"]
            nw_bus = network["bus"]["$(bus_id)"]

            if gen["qg"] < nw_gen["qmax"] && gen["qg"] > nw_gen["qmin"]
                bus = cont_sol["bus"]["$(bus_id)"]
                #if !isapprox(bus["vm"], nw_bus["vm"])
                if !isapprox(bus["vm"], nw_bus["vm"], atol=vm_eq_tol/2)
                    #debug(LOGGER, "bus $(bus_id) : vm_base $(nw_bus["vm"]) - vm $(bus["vm"]) : reactive bounds $(nw_gen["qmin"]) - $(gen["qg"]) - $(nw_gen["qmax"])")
                    warn(LOGGER, "update vm on bus $(bus_id) in contingency $(label) to match base case $(bus["vm"]) -> $(nw_bus["vm"]) due to within reactive bounds")
                end
                bus["vm"] = nw_bus["vm"]
            end

            pg_calc = nw_gen["pg"] + nw_gen["alpha"]*delta
            pg_calc = max(pg_calc, nw_gen["pmin"])
            pg_calc = min(pg_calc, nw_gen["pmax"])

            if !isapprox(gen["pg"], pg_calc, atol=1e-5)
                warn(LOGGER, "pg value on gen $(i) $(nw_gen["source_id"]) in contingency $(label) is not consistent with the computed value given:$(gen["pg"]) calc:$(pg_calc)")
            end

            if gen["pg"] > nw_gen["pmax"]
                warn(LOGGER, "update pg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match ub $(gen["pg"]) -> $(nw_gen["pmax"])")
                push!(pg_changes, gen["pg"] - nw_gen["pmax"])
                gen["pg"] = nw_gen["pmax"]
            end

            if gen["pg"] < nw_gen["pmin"]
                warn(LOGGER, "update pg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match lb $(gen["pg"]) -> $(nw_gen["pmin"])")
                push!(pg_changes, nw_gen["pmin"] - gen["pg"])
                gen["pg"] = nw_gen["pmin"]
            end

            if gen["qg"] > nw_gen["qmax"]
                warn(LOGGER, "update qg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match ub $(gen["qg"]) -> $(nw_gen["qmax"])")
                push!(qg_changes, gen["qg"] - nw_gen["qmax"])
                gen["qg"] = nw_gen["qmax"]
            end

            if gen["qg"] < nw_gen["qmin"]
                warn(LOGGER, "update qg on gen $(i) $(nw_gen["source_id"]) in contingency $(label) to match lb $(gen["qg"]) -> $(nw_gen["qmin"])")
                push!(qg_changes, nw_gen["qmin"] - gen["qg"])
                gen["qg"] = nw_gen["qmin"]
            end
        else
            gen["pg"] = 0.0
            gen["qg"] = 0.0
        end
    end

    # test imbalance on branch conts after flow corrections
    #=
    network_tmp = deepcopy(network)
    cont_type = cont_sol["cont_type"]
    cont_idx = cont_sol["cont_comp_id"]
    network_tmp["branch"]["$(cont_idx)"]["br_status"] = 0
    PowerModels.update_data!(network_tmp, cont_sol)
    deltas = compute_power_balance_deltas!(network_tmp)
    println("$(label): $(deltas)")
    =#

    cont_changed = length(vm_changes) > 1 || length(bs_changes) > 1 || length(pg_changes) > 1 || length(qg_changes) > 1

    if cont_changed
        _summary_changes(network, label, vm_changes, bs_changes, pg_changes, qg_changes)
    end

    return (changed=Int(cont_changed), vm_changes_max=maximum(vm_changes), bs_changes_max=maximum(bs_changes), pg_changes_max=maximum(pg_changes), qg_changes_max=maximum(qg_changes))
end


@everywhere function _summary_changes(network, contingency, vm_changes, bs_changes, pg_changes, qg_changes)
    println("")

    data = [
        "----",
        "contingency",
        "bus",
        "branch",
        "gen_cont",
        "branch_cont",
        "vm_count",
        "bs_count",
        "pg_count",
        "qg_count",
        "vm_max",
        "bs_max",
        "pg_max",
        "qg_max",
        "vm_mean",
        "bs_mean",
        "pg_mean",
        "qg_mean",
    ]
    println(join(data, ", "))

    data = [
        "DATA_CHANGES",
        contingency,
        length(network["bus"]),
        length(network["branch"]),
        length(network["gen_contingencies"]),
        length(network["branch_contingencies"]),
        length(vm_changes)-1,
        length(bs_changes)-1,
        length(pg_changes)-1,
        length(qg_changes)-1,
        maximum(vm_changes),
        maximum(bs_changes),
        maximum(pg_changes),
        maximum(qg_changes),
        mean(vm_changes),
        mean(bs_changes),
        mean(pg_changes),
        mean(qg_changes),
    ]
    println(join(data, ", "))
end


@everywhere function write_solution2_contingency(io::IO, pm_network, contingency_solution)
    base_mva = pm_network["baseMVA"]

    bus_switched_shunt_b = Dict(i => 0.0 for (i,bus) in pm_network["bus"])
    for (i,nw_shunt) in pm_network["shunt"]
        if nw_shunt["dispatchable"] && nw_shunt["status"] == 1
            #@assert nw_shunt["gs"] == 0.0
            shunt = contingency_solution["shunt"][i]
            bus_switched_shunt_b["$(nw_shunt["shunt_bus"])"] += shunt["bs"]
        end
    end

    write(io, "-- contingency\n")
    write(io, "label\n")
    write(io, "$(contingency_solution["label"])\n")

    write(io, "-- bus section\n")
    write(io, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
    for (i,bus) in contingency_solution["bus"]
        nw_bus = pm_network["bus"][i]
        write(io, "$(nw_bus["index"]), $(bus["vm"]), $(rad2deg(bus["va"])), $(base_mva*bus_switched_shunt_b[i])\n")
    end

    write(io, "-- generator section\n")
    write(io, "i, id, p(MW), q(MVAR)\n")
    for (i,gen) in contingency_solution["gen"]
        nw_gen = pm_network["gen"][i]
        bus_index = nw_gen["source_id"][2]
        gen_id = nw_gen["source_id"][3]
        write(io, "$(bus_index), $(gen_id), $(base_mva*gen["pg"]), $(base_mva*gen["qg"])\n")
    end

    write(io, "-- delta section\n")
    write(io, "delta(MW)\n")
    write(io, "$(base_mva*contingency_solution["delta"])\n")
end


"checks feasibility criteria of network solution, produces an error if a problem is found"
function check_network_solution(network)
    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            if bus["vm"] > bus["vmax"] || bus["vm"] < bus["vmin"]
                error(LOGGER, "vm on $(bus["source_id"]) is not in bounds $(bus["vmin"]) to $(bus["vmax"]), given $(bus["vm"])")
            end
        end
    end

    for (i,shunt) in network["shunt"]
        if shunt["status"] != 0
            if haskey(shunt, "dispatchable")
                if shunt["dispatchable"]
                    @assert shunt["gs"] == 0.0
                    @assert haskey(shunt, "bmin") && haskey(shunt, "bmax")
                    if shunt["bs"] > shunt["bmax"] || shunt["bs"] < shunt["bmin"]
                        error(LOGGER, "bs on $(shunt["source_id"]) is not in bounds $(shunt["bmin"]) to $(shunt["bmax"]), given $(shunt["bs"])")
                    end
                end
            end
        end
    end

    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            if gen["pg"] > gen["pmax"] || gen["pg"] < gen["pmin"]
                error(LOGGER, "pg on gen $(gen["source_id"]) not in bounds $(gen["pmin"]) to $(gen["pmax"]), given $(gen["pg"])")
            end

            if gen["qg"] > gen["qmax"] || gen["qg"] < gen["qmin"]
                error(LOGGER, "pg on gen $(gen["source_id"]) not in bounds $(gen["qmin"]) to $(gen["qmax"]), given $(gen["qg"])")
            end
        end
    end
end


function compute_power_balance_deltas!(network)
    flows = PowerModels.calc_branch_flow_ac(network)
    PowerModels.update_data!(network, flows)

    balance = PowerModels.calc_power_balance(network)
    PowerModels.update_data!(network, balance)

    p_delta_abs = [abs(bus["p_delta"]) for (i,bus) in network["bus"] if bus["bus_type"] != 4]
    q_delta_abs = [abs(bus["q_delta"]) for (i,bus) in network["bus"] if bus["bus_type"] != 4]

    return (
        p_delta_abs_max = maximum(p_delta_abs),
        p_delta_abs_mean = mean(p_delta_abs),
        q_delta_abs_max = maximum(q_delta_abs),
        q_delta_abs_mean = mean(q_delta_abs),
    )
end


function combine_files(files, output_file_name; output_dir="")
    if length(output_dir) > 0
        output_path = joinpath(output_dir, output_file_name)
    else
        output_path = output_file_name
    end

    open(output_path, "w") do output
        for file in files
            open(file, "r") do input
                for line in readlines(input, keep=true)
                    write(output, line)
                end
            end
        end
    end

    return output_path
end


function remove_files(files)
    for file in files
        if isfile(file)
            info(LOGGER, "deleting: $(file)")
            rm(file)
        else
            info(LOGGER, "skipping file: $(file)")
        end
    end
end
