# PowerModelsSecurityConstrained.jl Documentation

```@meta
CurrentModule = PowerModelsSecurityConstrained
```

## Overview

This package provides an extension to PowerModels for Security-Constrained
Optimization problems.  The original implementation developed an SCOPF solver
that was used as the benchmark algorithm for ARPA-e's grid optimization
competition challenge 1 in October 2019.


## Installation

The latest stable release of PowerModelsSecurityConstrained can be installed using the Julia package manager with

```julia
] add PowerModelsSecurityConstrained
```

For the current development version, "checkout" this package with

```julia
] add PowerModelsSecurityConstrained#master
```

At least one solver is required for running PowerModelsSecurityConstrained.  The open-source solver Ipopt is recommended, as it is fast, scaleable and can be used to solve a wide variety of the problems and network formulations provided in PowerModels.  The Ipopt solver can be installed via the package manager with

```julia
] add Ipopt
```

Test that the package works by running

```julia
] test PowerModelsSecurityConstrained
```
