using Documenter, PowerModelsSecurityConstrained

makedocs(
    modules = [PowerModelsSecurityConstrained],
    sitename = "PowerModelsSecurityConstrained",
    authors = "Carleton Coffrin",
    pages = [
        "Home" => "index.md",
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/PowerModelsSecurityConstrained.jl.git",
)
