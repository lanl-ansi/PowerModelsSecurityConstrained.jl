module PowerModelsSecurityConstrained

using Distributed
using SparseArrays
import Statistics: mean

using Memento
using JSON

using JuMP

using PowerModels
using InfrastructureModels


const LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(PowerModels)`
# NOTE: If this line is not included then the precompiled `PowerModels._LOGGER` won't be registered at runtime.
__init__() = Memento.register(LOGGER)


const vm_eq_tol = 1e-4
const vm_bound_tol = 1e-4
const qg_bound_tol = 1e-4
const pg_loss_tol = 1e-6


"Suppresses information and warning messages output by PMSC, for fine grained control use the Memento package"
function silence()
    Memento.info(LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
    Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
    Memento.setlevel!(Memento.getlogger(PowerModels), "error")
    Memento.setlevel!(Memento.getlogger(PowerModelsSecurityConstrained), "error")
end


include("core/common.jl")
include("core/variable.jl")
include("core/constraint.jl")
include("core/constraint_template.jl")

include("form/apo.jl")

include("io/parsers.jl")

include("prob/opf.jl")
include("prob/scopf.jl")
include("prob/security-stage.jl")

include("util/contingency-filters.jl")

# this must come last to support automated export
include("core/export.jl")

end