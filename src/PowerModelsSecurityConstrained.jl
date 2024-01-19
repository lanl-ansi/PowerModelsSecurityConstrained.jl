module PowerModelsSecurityConstrained

using Distributed
using SparseArrays
import Statistics: mean

using Memento
using JSON

import JuMP
import JuMP: @variable, @constraint, @objective, @expression

import InfrastructureModels
const _IM = InfrastructureModels

import PowerModels; const _PM = PowerModels
import PowerModels: ids, ref, var, con, sol, nw_ids, nw_id_default



const _LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level _LOGGER at runtime so that folks can access the
# _LOGGER via `get_LOGGER(PowerModelsSecurityConstrained)`
# NOTE: If this line is not included then the precompiled
# `PowerModelsSecurityConstrained.__LOGGER` won't be registered at runtime.
function __init__()
   Memento.register(_LOGGER)
   _LOGGER.name = "PMSC" # note must come after register, see discussion in issue #17
end


"Suppresses information and warning messages output by PMSC, for fine grained control use the Memento package"
function silence()
    Memento.info(_LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
    Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
    Memento.setlevel!(Memento.getlogger(PowerModels), "error")
    Memento.setlevel!(Memento.getlogger(PowerModelsSecurityConstrained), "error")
end


include("core/data.jl")
include("core/variable.jl")
include("core/constraint_template.jl")
include("core/constraint.jl")
include("core/expression_template.jl")
include("core/expression.jl")
include("core/objective.jl")
include("core/ref.jl")

include("form/acp.jl")
include("form/acr.jl")
include("form/wr.jl")
include("form/apo.jl")
include("form/dcp.jl")

include("io/goc.jl")
include("io/goc_c1.jl")
include("io/goc_c2.jl")

include("prob/opf.jl")
include("prob/ots.jl")
include("prob/scopf.jl")
include("prob/contingency-stage.jl")

include("util/contingency-filters.jl")
include("util/scopf-iterative.jl")

# this must come last to support automated export
include("core/export.jl")

end