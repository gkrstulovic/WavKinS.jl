push!(LOAD_PATH, "../")
using Documenter, WavKinS


makedocs(
        sitename="WavKinS.jl",
        modules=[WavKinS],
        pages=[
                "Home" => "index.md",
                "Quick start" => [
                        "Installing WavKinS" => "quick_start/installation.md",
                        "Computing a collisional integral" => "quick_start/collisional_integral_tutorial.md",
                        "Running an example" => "quick_start/running_tutorial.md",
                        "Adding a diagnostic" => "quick_start/diagnostics_tutorial.md",
                ],
                "Basics" => "basics.md",
                "Physical systems" => [
                        "Acoustic" => "physical_systems/Acoustic.md",
                        "Bogoliubov" => "physical_systems/Bogoliubov.md",
                        "NLS" => "physical_systems/NLS.md",
                        "MMT" => "physical_systems/MMT.md",
                        "Petviashvilli" => "physical_systems/Petviashvilli.md",
                        "Smoluchowski" => "physical_systems/Smoluchowski.md",
                        "Stratified" => "physical_systems/Stratified.md",
                        "Generalities" => "physical_systems/generalities.md",
                        "TestTimeStepping" => "physical_systems/TestTimeStepping.md"
                ],
                "Numerical methods" => "numerical-methods.md",
                "Multithreading" => "multithreading.md",
                "Plots" => "plots.md",
                "Creating new solvers" => "creating-new-solvers.md"
        ],
        format=Documenter.HTML(prettyurls=false),
        remotes=nothing
)
