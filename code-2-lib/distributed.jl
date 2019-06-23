using Distributed

# this appears to be the point of diminishing returns
max_node_processes = 12

"turns a slurm node list into a list of node names"
function parse_node_names(slurm_node_list::String)
    node_names = String[]
    slurm_node_list_parts = split(slurm_node_list, "[")
    node_name_prefix = slurm_node_list_parts[1]

    #if only one node was found then this length will be 1
    if length(slurm_node_list_parts) == 1
        push!(node_names, node_name_prefix)
    else
        # this assumes all nodes have a common prefix and the nodelist has the form [1-2,4,5-10,12]
        node_ids = split(split(slurm_node_list_parts[2], "]")[1], ",")
        for node_id in node_ids
            node_id_parts = split(node_id, "-")
            start = parse(Int, node_id_parts[1])
            stop = parse(Int, length(node_id_parts) == 1 ? node_id_parts[1] : node_id_parts[2])
            digits = length(node_id_parts[1])
            for i in start:stop
                node_name = "$(node_name_prefix)$(lpad(string(i), digits, "0"))"
                push!(node_names, node_name)
            end
        end
    end
    return node_names
end

function add_local_procs()
    proc_ids = Int[]
    if Distributed.nprocs() <= 2
        node_processes = min(trunc(Int, Sys.CPU_THREADS*0.75), max_node_processes)
        println("local processes: $(node_processes) of $(Sys.CPU_THREADS)")

        proc_ids = Distributed.addprocs(node_processes)
    end
    return proc_ids
end

function add_remote_procs()
    if haskey(ENV, "SLURM_JOB_NODELIST")
        node_list = ENV["SLURM_JOB_NODELIST"]
    elseif haskey(ENV, "SLURM_NODELIST")
        node_list = ENV["SLURM_NODELIST"]
    else
        @info("unable to find slurm node list environment variable")
        return Int[]
    end
    #println(node_list)

    node_names = parse_node_names(node_list)
    #println(node_names)

    node_names = [name for name in node_names if name != gethostname()]

    proc_ids = Int[]
    if length(node_names) > 0
        node_processes = min(trunc(Int, Sys.CPU_THREADS*0.75), max_node_processes)
        println("remote processes per node: $(node_processes)/$(Sys.CPU_THREADS)")
        for i in 1:node_processes
            node_proc_ids = Distributed.addprocs(node_names, sshflags="-oStrictHostKeyChecking=no")
            println("process id batch $(i) of $(node_processes): $(node_proc_ids)")
            for npid in node_proc_ids
                push!(proc_ids, npid)
            end
        end
    else
        @info("no remote slurm nodes found")
    end

    #println("process ids: $(proc_ids)")
    return proc_ids
end


function add_procs()
    proc_ids = add_local_procs()
    remote_proc_ids = add_remote_procs()

    for pid in remote_proc_ids
        push!(proc_ids, pid)
    end

    println("process ids: $(proc_ids)")
    return proc_ids
end
