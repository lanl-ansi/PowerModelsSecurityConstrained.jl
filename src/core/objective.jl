
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

