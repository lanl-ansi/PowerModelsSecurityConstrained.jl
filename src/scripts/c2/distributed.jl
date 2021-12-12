using Distributed

# this appears to be the point of diminishing returns on most harware platforms
max_node_processes = 32

node_processes = min(Sys.CPU_THREADS - 2, max_node_processes)

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
        @info("remote processes per node: $(node_processes)/$(Sys.CPU_THREADS)")
        tasks = []
        for i in 1:node_processes
            task = @async Distributed.addprocs(node_names, topology=:master_worker, sshflags="-oStrictHostKeyChecking=no")
            push!(tasks, task)
            #node_proc_ids = Distributed.addprocs(node_names, topology=:master_worker, sshflags="-oStrictHostKeyChecking=no")
            #@info("process id batch $(i) of $(node_processes): $(node_proc_ids)")
        end
    for (i,task) in enumerate(tasks)
            wait(task)
            @info("task $(i) of $(node_processes) complete")
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




function check_available_memory(bytes_required; memory_margin::Real=2.5)
    start = time()

    if haskey(ENV, "SLURM_JOB_NODELIST")
        node_list = ENV["SLURM_JOB_NODELIST"]
    elseif haskey(ENV, "SLURM_NODELIST")
        node_list = ENV["SLURM_NODELIST"]
    else
        node_list = "none"
    end

    node_count = length(parse_node_names(node_list))

    workers = Distributed.workers()

    required_mem = memory_margin*length(workers)*bytes_required

    if required_mem > node_count*Sys.total_memory()
        gb_scale = 2^30
        task_mem_gb = trunc(memory_margin*bytes_required/gb_scale, digits=2)
        required_mem_gb = trunc(required_mem/gb_scale/node_count, digits=2)
        system_mem_gb = trunc(Sys.total_memory()/gb_scale, digits=2)
        @info("insufficient memory to run $(length(workers)) workers with $(task_mem_gb) Gb per worker per node; required $(required_mem_gb) Gb per node, node system $(system_mem_gb) Gb")

        worker_max = node_count*floor(Int, Sys.total_memory()/(memory_margin*bytes_required))
        workers_to_remove = length(workers) - worker_max
        @info("removing $(workers_to_remove) workers across $(node_count) nodes")

        step_size = length(workers)/workers_to_remove
        idx = 1.0
        while idx < length(workers)
            worker_id = workers[floor(Int, idx)]
            rmprocs(worker_id)
            idx += step_size
        end

        @info("worker ids: $(Distributed.workers())")
    end

    @info("worker revise time: $(time() - start)")
end
