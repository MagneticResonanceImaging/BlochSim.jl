tests = [
    "helpers",
    "Spin",
    "sequences"
]
for t in tests
  include("$(t).jl")
end
