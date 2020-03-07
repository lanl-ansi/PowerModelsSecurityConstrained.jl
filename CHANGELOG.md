PowerModelsSecurityConstrained.jl Change Log
============================================

### Staged
- Generalized dc-only OPF formulations
- Updated `bsh` variable name to `bs`

### Staged
- Added variant of the `build_scopf_dc_cuts_soft` formulation

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
