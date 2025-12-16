push!(LOAD_PATH, "../src/")
using Documenter, Dynamic_objectives

makedocs(
    sitename = "Dynamic_objectives Documentation",
    modules = [Dynamic_objectives],
    format = Documenter.HTML(
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Architecture" => "architecture.md",
        "Models" => [
            "Model Catalog" => "model_catalog.md",
            "Testing Guide" => "testing_guide.md"
        ],
        "Integration" => [
            "Globtim Integration" => "integration.md",
            "PostProcessing Requirements" => "postprocessing_requirements.md"
        ],
        "API Reference" => "api_reference.md",
        "Performance Report" => "performance_report.md"
    ],
    checkdocs = :none
)
