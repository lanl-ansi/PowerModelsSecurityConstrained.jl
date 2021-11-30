# PMSC Scripts

These scripts implement the SCOPF solver that was used in the final event of
the grid optimization competition challenge 1 in October 2019.  The primary
entry point for the solver is `goc_challenge1_huristic_cli.jl`.

These scripts also include the tooling to run the solver on the optimization
competition platform, i.e. `MyJulia1.jl` and `MyJulia2.jl`.

Note that `goc_challenge1_huristic_cli.jl` is configured to use`Cbc` for
solving LP problems.  Commercial solvers like Gurobi and CPlex are highly
recommended if a license is available.  Gurobi v8 was used by this software
during the final event of the grid optimization competition challenge 1.

