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
    if Distributed.nprocs() >= 2
        return
    end
    node_processes = min(trunc(Int, Sys.CPU_THREADS*0.75), max_node_processes)
    @info("local processes: $(node_processes) of $(Sys.CPU_THREADS)")

    Distributed.addprocs(node_processes, topology=:master_worker)
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

    @info("host name: $(gethostname())")
    @info("slurm allocation nodes: $(node_names)")

    # some systems add .local to the host name
    hostname = gethostname()
    if endswith(hostname, ".local")
        hostname = hostname[1:end-6]
    end
    if endswith(hostname, ".localdomain")
        hostname = hostname[1:end-12]
    end

    node_names = [name for name in node_names if name != hostname]

    if length(node_names) > 0
        @info("remote slurm nodes: $(node_names)")
        node_processes = min(trunc(Int, Sys.CPU_THREADS*0.75), max_node_processes)
        @info("remote processes per node: $(node_processes)/$(Sys.CPU_THREADS)")
        for i in 1:node_processes
            node_proc_ids = Distributed.addprocs(node_names, topology=:master_worker, sshflags="-oStrictHostKeyChecking=no")
            @info("process id batch $(i) of $(node_processes): $(node_proc_ids)")
        end
    else
        @info("no remote slurm nodes found")
    end
end


function add_procs()
    start = time()
    add_local_procs()
    @info("local proc start: $(time() - start)")

    start = time()
    add_remote_procs()
    @info("remote proc start: $(time() - start)")

    @info("worker ids: $(Distributed.workers())")
end
