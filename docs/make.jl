using Documenter
using DynamicObjectives

makedocs(
    sitename = "DynamicObjectives.jl",
    format = Documenter.HTML(
        prettyurls = false,
        canonical = "https://gescholt.github.io/DynamicObjectives.jl",
        # The API page is one @autodocs block over ~173 docstrings, which lands
        # around 270 KiB -- over Documenter's 200 KiB default.
        size_threshold = 400 * 1024,
        size_threshold_warn = 300 * 1024,
    ),
    pages = [
        "Home" => "index.md",
        "Model Catalog" => "model_catalog.md",
        "Glued Objectives" => "glued_objectives.md",
        "API Reference" => "api.md",
    ],
    modules = [DynamicObjectives],
)

deploydocs(
    repo = "github.com/gescholt/DynamicObjectives.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
)
