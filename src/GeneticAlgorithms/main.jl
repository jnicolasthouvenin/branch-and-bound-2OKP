
# This module contains the implementation of two EMOA algorithms: NSGAII and SMS-EMOA
# Implementation follows the general algorithms presented in:
# - source nsgaii  : https://doi.org/10.1007/3-540-45356-3_83
# - source sms-emoa: https://doi.org/10.1007/978-3-540-31880-4_5

# This project was part of the course Multi-objective Metaheuristics from the master 2 ORO at University of Nantes

module GeneticAlgorithms

# Libraries
using Random

# Structs, tools, debuging, indicators computation...
include("structs.jl")
include("printing.jl")
include("debug.jl")
include("utils.jl")
include("solution.jl")
include("indicator.jl")

# GA operators
include("initialisation.jl")
include("selection.jl")
include("crossover.jl")
include("mutation.jl")
include("reparation.jl")

# EMOA algorithms
include("nsgaii.jl")
include("sms_emoa.jl")

# INTRO #####################################################

debug = false       # wether of not assertions are checked (expensive, use only for debugging)

# CODE #######################################################

# returns the configuration considered best by our experiments
function get_optimal_config(prob::_MOMKP, algorithm)
    size, crossover, mutation_rate = 0,0,0

    if algorithm == NSGAII_solve || algorithm == SMS_EMOA_solve
        if prob.n < 50
            size = 75
            crossover = CROSSOVER_binary_uniform
            mutation_rate = 0.05
        elseif prob.n < 150
            size = 150
            crossover = CROSSOVER_binary_uniform
            mutation_rate = 0.05
        else
            size = 150
            crossover = CROSSOVER_two_points
            mutation_rate = 0.2
        end
    else
        println("ERROR : algorithm not recognised")
        return false
    end

    config = Config(
        size,
        crossover,
        mutation_rate,
        REPARATION_advanced!
    )

    return config
end

function solve(prob::_MOMKP; max_t = 15.0, algorithm_solve = NSGAII_solve, seeding_sols = Vector{Sol}(undef, 0))

    config = get_optimal_config(prob, algorithm_solve)
    front = algorithm_solve(prob, config, seeding_sols = seeding_sols, max_t = max_t, verbose = false)

    return front

end

end