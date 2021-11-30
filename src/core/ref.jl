
""
function ref_c1!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    _PM.apply_pm!(_ref_c1!, ref, data)
end

function _ref_c1!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
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
