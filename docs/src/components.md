# Solver Components

```@meta
CurrentModule = PowerModelsSecurityConstrained
```

## Data Wangling

These tools are used for parsing data files and translating data into the
PowerModels internal data model.

```@docs
parse_goc_files
build_pm_model
```

## OPF Formulations

These are standard OPF formulations (i.e. without contingency constraints) and
are used as sub-processes in building contingency constrained solutions.

```@docs
run_opf_shunt
run_opf_cheap
run_opf_cheap_lazy_acr
```

## Contingency Filters

These are tools for checking for constraint violations in contingencies

```@docs
check_contingencies_branch_power
check_contingencies_branch_power_bpv
```

## SCOPF Formulations

These are contingency constrained OPF formulations.

```@docs
run_scopf_cuts_soft
run_scopf_cuts_soft_bpv
```

## Second-Stage Solvers

These are solvers for quickly solving large collections of contingencies.
These solvers are usually more detailed on than the models used in the
contingency filters.

```@docs
run_fixpoint_pf_pvpq!
run_fixpoint_pf_bqv!
```
