PowerModelsSecurityConstrained.jl Change Log
============================================

### Staged
- nothing

### v0.4.0
- Added support for Memento v0.13, v1.0
- Added `expression_branch_power_yt_from/to` for common branch flow expressions (#2)
- Improved generality of contingency problem definitions (#2)
- Improved power balance constraint resilience
- Generalized dc-only OPF and SCOPF formulations (#2)
- Updated `bsh` variable name to `bs` (#5)
- Updated scope of data processing functions (#4)
- Updated logger printed name to PMSC (#17)
- Reorganize code to standard PowerModels locations (#4)
- Removed duplicate `objective_variable_pg_cost` implementation
- Remove PowerModel parameters from fixed-formulation models
- Renamed problem functions to have less cryptic names (PR #18)

### v0.3.2
- Add variant of the `build_scopf_dc_cuts_soft` formulation

### v0.3.1
- Added support for JuMP v0.21

### v0.3.0
- Updated to PowerModels v0.15
- Changed post_* function names to build_* to follow new naming conventions

### v0.2.1
- Add support for PowerModels v0.14

### v0.2.0
- Updated to PowerModels v0.13
- Removed try/catch blocks around solver calls
- Removed duplicate psse data parser (#3)

### v0.1.1
- Added missing `gen_default` dictionary for security-stage solvers
- Added tests for infeasible network models
- Improved the resilience of `goc_challenge1_huristic` to infeasible models
- Fixed contingency area indexing bug

### v0.1.0
- ARPA-e Grid Optimization Competition Challenge 1, Final Event submission
