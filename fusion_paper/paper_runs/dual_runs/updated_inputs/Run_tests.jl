cases_2_run = [
    "test_updated_basecase_cost8500/Run.jl",
    "test_updated_basecase_cost3000/Run.jl",
    "test_updated_basecase_cost8500_thermstor/Run.jl",
]

for case in cases_2_run
    println("Running case: ", case)
    include(joinpath(@__DIR__, case))
end