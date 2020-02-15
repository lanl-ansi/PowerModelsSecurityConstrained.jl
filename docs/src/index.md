# PowerModelsSecurityConstrained.jl Documentation

```@meta
CurrentModule = PowerModelsSecurityConstrained
```

## Overview


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
