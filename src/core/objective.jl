"""
Checks that all generator cost models are of the same type

adapted from the implementation in PowerModels <= v0.19
"""
function _check_gen_cost_models(pm::_PM.AbstractPowerModel)
    model = nothing

    for (n, nw_ref) in _PM.nws(pm)
        for (i,gen) in nw_ref[:gen]
            if haskey(gen, "cost")
                if model == nothing
                    model = gen["model"]
                else
                    if gen["model"] != model
                        Memento.error(_LOGGER, "cost models are inconsistent, the typical model is $(model) however model $(gen["model"]) is given on generator $(i)")
                    end
                end
            else
                Memento.error(_LOGGER, "no cost given for generator $(i)")
            end
        end
    end

    return model
end


"""
adds pg_cost variables and constraints

adapted from the implementation in PowerModels <= v0.19
"""
function _objective_variable_pg_cost(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            pg_var = var(pm, n, :pg, i)
            pmin = JuMP.lower_bound(pg_var)
            pmax = JuMP.upper_bound(pg_var)

            points = _PM.calc_pwl_points(gen["ncost"], gen["cost"], pmin, pmax)

            pg_cost_lambda = JuMP.@variable(pm.model,
                [i in 1:length(points)], base_name="$(n)_pg_cost_lambda",
                lower_bound = 0.0,
                upper_bound = 1.0
            )
            JuMP.@constraint(pm.model, sum(pg_cost_lambda) == 1.0)

            pg_expr = 0.0
            pg_cost_expr = 0.0
            for (i,point) in enumerate(points)
                pg_expr += point.mw*pg_cost_lambda[i]
                pg_cost_expr += point.cost*pg_cost_lambda[i]
            end
            JuMP.@constraint(pm.model, pg_expr == pg_var)
            pg_cost[i] = pg_cost_expr
        end

        report && _PM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end


""
function objective_c1_variable_pg_cost(pm::_PM.AbstractPowerModel; kwargs...)
    model = _check_gen_cost_models(pm)

    if model == 1
        return _objective_variable_pg_cost(pm; kwargs...)
    elseif model == 2
        return objective_variable_pg_cost_polynomial_linquad(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_polynomial_linquad(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (nw, nw_ref) in _PM.nws(pm)
        pg_cost = var(pm, nw)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, nw, :gen)
            pg = var(pm, nw, :pg, i)

            if length(gen["cost"]) == 1
                pg_cost[i] = gen["cost"][1]
            elseif length(gen["cost"]) == 2
                pg_cost[i] = gen["cost"][1]*pg + gen["cost"][2]
            elseif length(gen["cost"]) == 3
                pg_cost[i] = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
            else
                pg_cost[i] = 0.0
            end
        end

        report && _PM.sol_component_value(pm, nw, :gen, :pg_cost, ids(pm, nw, :gen), pg_cost)
    end
end



""
function objective_c1_variable_pg_cost_basecase(pm::_PM.AbstractPowerModel; kwargs...)
    model = _check_gen_cost_models(pm)

    if model == 1
        return objective_c1_variable_pg_cost_basecase_pwl(pm; kwargs...)
    elseif model == 2
        return objective_c1_variable_pg_cost_basecase_polynomial_linquad(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

"adds pg_cost variables and constraints for base case only"
function objective_c1_variable_pg_cost_basecase_pwl(pm::_PM.AbstractPowerModel, nw::Int=0, report::Bool=true)
    pg_cost = var(pm, nw)[:pg_cost] = Dict{Int,Any}()

    for (i,gen) in ref(pm, nw, :gen)
        pg_var = var(pm, nw, :pg, i)
        pmin = JuMP.lower_bound(pg_var)
        pmax = JuMP.upper_bound(pg_var)

        # note pmin/pmax may be different from gen["pmin"]/gen["pmax"] in the on/off case
        points = _PM.calc_pwl_points(gen["ncost"], gen["cost"], pmin, pmax)

        pg_cost_lambda = JuMP.@variable(pm.model,
            [i in 1:length(points)], base_name="$(nw)_pg_cost_lambda",
            lower_bound = 0.0,
            upper_bound = 1.0
        )
        JuMP.@constraint(pm.model, sum(pg_cost_lambda) == 1.0)

        pg_expr = 0.0
        pg_cost_expr = 0.0
        for (i,point) in enumerate(points)
            pg_expr += point.mw*pg_cost_lambda[i]
            pg_cost_expr += point.cost*pg_cost_lambda[i]
        end
        JuMP.@constraint(pm.model, pg_expr == pg_var)
        pg_cost[i] = pg_cost_expr
    end

    report && _PM.sol_component_value(pm, nw, :gen, :pg_cost, ids(pm, nw, :gen), pg_cost)
end


"adds pg_cost variables and constraints for base case only"
function objective_c1_variable_pg_cost_basecase_polynomial_linquad(pm::_PM.AbstractPowerModel, nw::Int=0, report::Bool=true)
    pg_cost = var(pm, nw)[:pg_cost] = Dict{Int,Any}()

    for (i,gen) in ref(pm, nw, :gen)
        pg = var(pm, nw, :pg, i)

        if length(gen["cost"]) == 1
            pg_cost[i] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            pg_cost[i] = gen["cost"][1]*pg + gen["cost"][2]
        elseif length(gen["cost"]) == 3
            pg_cost[i] = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
        else
            pg_cost[i] = 0.0
        end
    end

    report && _PM.sol_component_value(pm, nw, :gen, :pg_cost, ids(pm, nw, :gen), pg_cost)
end


"adds pd_value variables and constraints"
function objective_c2_variable_pd_value(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        pd_value = var(pm, n)[:pd_value] = Dict{Int,Any}()

        for (i,load) in ref(pm, n, :load)
            pd = load["pd_nominal"]*var(pm, n, :z_demand, i)

            if load["model"] == 1
                points = _PM.calc_pwl_points(load["ncost"], load["cost"], load["pd_min"], load["pd_max"])

                if length(points) == 2 && isapprox(points[1].cost, points[2].cost)
                    #println("c")
                    pd_value[i] = points[1].cost
                else
                    pd_value_lambda = JuMP.@variable(pm.model,
                        [i in 1:length(points)], base_name="$(n)_pd_value_lambda",
                        lower_bound = 0.0,
                        upper_bound = 1.0
                    )
                    JuMP.@constraint(pm.model, sum(pd_value_lambda) == 1.0)

                    pd_expr = 0.0
                    pd_value_expr = 0.0
                    for (i,point) in enumerate(points)
                        pd_expr += point.mw*pd_value_lambda[i]
                        pd_value_expr += point.cost*pd_value_lambda[i]
                    end
                    JuMP.@constraint(pm.model, pd_expr == pd)
                    pd_value[i] = pd_value_expr
                end

            elseif load["model"] == 2
                if length(load["cost"]) == 1
                    pd_value[i] = load["cost"][1]
                elseif length(load["cost"]) == 2
                    pd_value[i] = load["cost"][1]*pd + load["cost"][2]
                elseif length(load["cost"]) == 3
                    pd_value[i] = load["cost"][1]*pd^2 + load["cost"][2]*pd + load["cost"][3]
                else
                    pd_value[i] = 0.0
                end

            else # cost model not defined
                @assert(false)
            end
        end

        report && _PM.sol_component_value(pm, n, :load, :pd_value, ids(pm, n, :load), pd_value)
    end
end

