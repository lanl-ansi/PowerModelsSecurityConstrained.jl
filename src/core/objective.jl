
"adds pg_cost variables and constraints for base case only"
function objective_c1_variable_pg_cost_basecase(pm::_PM.AbstractPowerModel, nw::Int=0, report::Bool=true)
    pg_cost = var(pm, nw)[:pg_cost] = Dict{Int,Any}()

    for (i,gen) in ref(pm, nw, :gen)
        pg_vars = [var(pm, nw, :pg, i)[c] for c in _PM.conductor_ids(pm, nw)]
        pmin = sum(JuMP.lower_bound.(pg_vars))
        pmax = sum(JuMP.upper_bound.(pg_vars))

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
        JuMP.@constraint(pm.model, pg_expr == sum(pg_vars))
        pg_cost[i] = pg_cost_expr
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

