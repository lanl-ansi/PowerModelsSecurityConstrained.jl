using Documenter, PowerModelsSecurityConstrained

makedocs(
    warnonly = Documenter.except(:linkcheck),
    modules = [PowerModelsSecurityConstrained],
    sitename = "PowerModelsSecurityConstrained",
    authors = "Carleton Coffrin",
    pages = [
        "Home" => "index.md",
        "Solver Components" => "components.md",
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/PowerModelsSecurityConstrained.jl.git",
)
