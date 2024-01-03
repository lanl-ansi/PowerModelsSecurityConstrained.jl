PowerModelsSecurityConstrained.jl Change Log
============================================

### Staged
- nothing

### v0.11.0
- Update to PowerModels v0.20
- Drop support for JuMP v0.22 and v0.23 (breaking)
- Drop support for JSON v0.18, v0.19, v0.20 (breaking)

### v0.10.1
- Add support for polynomial generator costs
- Fixed bug when linear solve failed during contingency checking
- Fixed PowerModels deprecation warnings (update `run_*` to `solve_*`)

### v0.10.0
- Add support for JuMP v1.0
- Drop support for JuMP v0.21
- Replace CBC with HiGHS in tests

### v0.9.1
- Add support for JuMP v0.23
- Update minimum Julia version to v1.6 (LTS)

### v0.9.0
- ARPA-e Grid Optimization Competition Challenge 2, Final Event submission (#34)
- Reorganize function names and file locations to support multiple competition solvers

### v0.8.2
- Add support for Memento v1.3

### v0.8.1
- Update to PowerModels v0.19

### v0.8.0
- Update to PowerModels v0.18

### v0.7.1
- Fix bug in `contingency_order` logic

### v0.7.0
- Added flexible multi-network SCOPF formulation (#19)
- Added `run_opf_cheap_target_acp` OPF formulation
- Updated contigency filters to use PowerModels PTDF cut tools (#21)
- Updated violation computations to support NaN reactive power values
- Improved distributed computation resilience to cluster configurations
- Simplified `compute_violations` functions (breaking)
- Fixed SCOPF formulations to allow shunt optimization

### v0.6.0
- Update to new function name convention of PowerModels v0.17 (breaking)
- Updated to PowerModels v0.17 (breaking)
- Added support for Memento v1.1

### v0.5.0
- Updated to PowerModels v0.16

### v0.4.1
- Updated InfrastructureModels from using to import
- Added use Gurobi option to cli interface
- Minor fixes to cli interface in scripts

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
- Removed PowerModel parameters from fixed-formulation models
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
