# PowerModelsSecurityConstrained Documentation

```@meta
CurrentModule = PowerModelsSecurityConstrained
```

## Overview

This package provides an extension to PowerModels for Security-Constrained
Optimization problems.  The core routines provided in this software have been used in the benchmark algorithm of ARPA-e's grid optimization competitions.


## Installation

The latest stable release of PowerModelsSecurityConstrained can be installed using the Julia package manager with

```julia
] add PowerModelsSecurityConstrained
```

For the current development version, "checkout" this package with

```julia
] add PowerModelsSecurityConstrained#master
```

At least one solver is required for running PowerModelsSecurityConstrained.  The open-source solver Ipopt is recommended, as it is fast, scalable and can be used to solve a wide variety of the problems and network formulations provided in PowerModels.  The Ipopt solver can be installed via the package manager with

```julia
] add Ipopt
```

Test that the package works by running

```julia
] test PowerModelsSecurityConstrained
```

## ARPA-e SCOPF Benchmark

PowerModelsSecurityConstrained includes a variety of tools for solving problems
with contingency constraints.  The `src/script` directory includes the
specific algorithms that were submitted to the final event in ARPA-e's Grid
Optimization Competition Challenge 1 in October 2019 and Challenge 2 in September 2021.
