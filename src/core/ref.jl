
function ref_add_goc!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    if _IM.ismultinetwork(data)
        nws_data = data["nw"]
    else
        nws_data = Dict("0" => data)
    end

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        ref = ref[:nw][nw_id]

        ref[:branch_sm_active] = Dict(x for x in ref[:branch] if haskey(x.second, "rate_a"))

        ref[:shunt_const] = Dict(x for x in ref[:shunt] if (!haskey(x.second, "dispatchable") || !x.second["dispatchable"]))
        ref[:shunt_var] = Dict(x for x in ref[:shunt] if (haskey(x.second, "dispatchable") && x.second["dispatchable"]))

        bus_shunts_const = Dict((i, []) for (i,bus) in ref[:bus])
        for (i,shunt) in ref[:shunt_const]
            push!(bus_shunts_const[shunt["shunt_bus"]], i)
        end
        ref[:bus_shunts_const] = bus_shunts_const

        bus_shunts_var = Dict((i, []) for (i,bus) in ref[:bus])
        for (i,shunt) in ref[:shunt_var]
            push!(bus_shunts_var[shunt["shunt_bus"]], i)
        end
        ref[:bus_shunts_var] = bus_shunts_var
    end
end
