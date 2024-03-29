# Solver Components

```@meta
CurrentModule = PowerModelsSecurityConstrained
```

## Data Wangling

These tools are used for parsing data files and translating data into the
PowerModels internal data model.

```@docs
parse_c1_case
build_c1_pm_model
parse_c2_case
build_c2_pm_model
```

## OPF Formulations

These are standard OPF formulations (i.e. without contingency constraints) and
are used as sub-processes in building contingency constrained solutions.

```@docs
run_c1_opf_shunt
run_c1_opf_cheap
run_c1_opf_cheap_lazy_acr
build_c2_opf_soft
build_c2_opf_soft_ctg
```

## OPF with UC Formulations

These are standard OPF formulations (i.e. without contingency constraints) with discrete variables for supporting the commitment of generation units.

```@docs
build_c2_opf_uc
```

## OTS Formulations

These are standard transmission switching formulations (i.e. without contingency constraints) and are used to find economic improving topology changes in a network.

```@docs
build_c2_ots_soft
```

## Contingency Filters

These are tools for checking for constraint violations in contingencies

```@docs
check_c1_contingency_violations
check_c1_contingencies_branch_power
check_c1_contingencies_branch_power_bpv
```

## SCOPF Formulations

These are contingency constrained OPF formulations.

```@docs
run_c1_scopf
run_c1_scopf_cuts_soft
run_c1_scopf_cuts_soft_bpv
run_c1_scopf_contigency_cuts
run_c1_scopf_ptdf_cuts!
```

## Contingency-Stage Solvers

These are solvers for quickly solving large collections of contingencies.
These solvers are usually more detailed on than the models used in the
contingency filters.

```@docs
run_c1_fixpoint_pf_pvpq!
run_c1_fixpoint_pf_bqv!
```
