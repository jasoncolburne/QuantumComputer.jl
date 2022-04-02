push!(LOAD_PATH,"../src/")

using Documenter, DocumenterMarkdown, Example

makedocs(
  sitename = "QuantumComputer.jl",
  format = Markdown()
)