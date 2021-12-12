# PMSC Scripts

These scripts implement the SCOPF benchmark solvers that were used in the final event of ARPA-e's grid optimization competition challenges. The directory `c1` includes the solver from Challenge 1 (October, 2019) and `c2` include the solver from Challenge 2 (September, 2021).  The primary entry point for each solver is `c1/goc_c1_huristic_cli.jl` and `c2/goc_c2_huristic_cli.jl` respectivly.

These scripts also include the tooling to run the solver on the grid optimization
competition's platform, i.e. `MyJulia1.jl` and `MyJulia2.jl`.

Note that `goc_c1_huristic_cli.jl` is configured to use `Cbc` for
solving LP problems. Commercial solvers like Gurobi and CPlex are highly
recommended if a license is available. Gurobi v8 was used by this software
during the final event of the grid optimization competition challenge 1.
