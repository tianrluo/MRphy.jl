using Documenter, MRphy

makedocs(
    modules = [MRphy],
    sitename="MRphy.jl",
    authors = "Tianrui Luo",
)

deploydocs(
    repo = "github.com/tianrluo/MRphy.jl.git"
)
