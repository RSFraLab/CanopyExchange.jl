using Documenter
using CanopyExchange

@info "Building docs with CanopyExchange source" path = pathof(CanopyExchange)

function build()
    pages = Any[
        "Home" => "index.md",
        "Concepts & Host Coupling" => [
            "Measures & Conservation" => "pages/measures_and_conservation.md",
            "CliMA Interoperability" => "pages/clima_interoperability.md",
            "Validation & Bit Parity" => "pages/validation.md",
        ],
        "API Reference" => [
            "Overview & Index" => "api/index.md",
            "Quantities & Layer Geometries" => "api/quantities_and_geometries.md",
            "Area Sources & Divergence" => "api/sources_and_tendencies.md",
            "CliMA Contracts & Scalar Traits" => "api/clima_and_scalars.md",
        ],
    ]

    mathengine = MathJax(
        Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict(),
            ),
        ),
    )

    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = nothing,
        repolink = nothing,
        mathengine = mathengine,
        collapselevel = 1,
        size_threshold = nothing,
        size_threshold_warn = nothing,
        example_size_threshold = nothing,
    )

    DocMeta.setdocmeta!(
        CanopyExchange,
        :DocTestSetup,
        :(using CanopyExchange);
        recursive = true,
    )

    makedocs(
        sitename = "CanopyExchange.jl",
        format = format,
        remotes = nothing,
        clean = true,
        doctest = true,
        checkdocs = :all,
        modules = [CanopyExchange],
        pages = pages,
    )
end

build()

if get(ENV, "CI", "false") == "true" && haskey(ENV, "GITHUB_REPOSITORY")
    deploydocs(
        repo = "github.com/$(ENV["GITHUB_REPOSITORY"]).git",
        target = "build",
        devbranch = "main",
        push_preview = true,
    )
end
