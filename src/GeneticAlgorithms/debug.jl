
# DEBUG SEVERAL SOLUTIONS #####################################################

function DEBUG_correct_solutions(solutions::Vector{Sol})
    for sol in solutions
        DEBUG_correct_solution(sol)
    end
end

# DEBUG ONE MOMKP SOLUTION

function DEBUG_correct_solution(sol::Sol)

    momkp = sol.prob
    n = momkp.n

    @assert length(sol.x) == n "sol.x not of the right size"

    z_true = [0.0, 0.0]
    w_true = zeros(Float64, momkp.m)
    for var in 1:n
        @assert sol.x[var] == 0 || sol.x[var] == 1 "sol.x has values different from 0 or 1, so not a knapsack solution."
        if sol.x[var] == 1
            z_true[1] += momkp.P[1,var]
            z_true[2] += momkp.P[2,var]
            for dim in 1:momkp.m
                w_true[dim] += momkp.W[dim,var]
            end
        end
    end

    @assert z_true[1] == sol.z[1] && z_true[2] == sol.z[2] "sol.z != z_true, $(sol.z) != $z_true"
    @assert w_true == sol.w "sol.w != w_true, $(sol.w) != $w_true"
end

# DEBUG SEVERAL FEASIBLE SOLUTIONS #####################################################

function DEBUG_feasible_solutions(solutions::Vector{Sol})
    for sol in solutions
        DEBUG_feasible_solution(sol)
    end
end

# DEBUG ONE FEASIBLE MOMKP SOLUTION

function DEBUG_feasible_solution(sol::Sol)
    DEBUG_correct_solution(sol)

    for dim in 1:sol.prob.m
        @assert sol.w[dim] <= sol.prob.ω[dim] "Solution not feasible, sol.w[$dim] > sol.prob.ω[$dim] : $(sol.w[dim]) > $(sol.prob.ω[dim])"
    end
end

# VARIOUS DEBUG TESTS #####################################################

# check if the solutions are mutually non dominated
function DEBUG_is_front(front::Vector{Sol})
    n = length(front)

    for i in 1:n
        for j in i+1:n
            # if one dominates the other, error
            @assert !((front[i].z[1] >= front[j].z[1] && front[i].z[2] > front[j].z[2]) || (front[i].z[1] > front[j].z[1] && front[i].z[2] >= front[j].z[2])) "All solutions in front aren't mutually non-dominated"
            @assert !((front[j].z[1] >= front[i].z[1] && front[j].z[2] > front[i].z[2]) || (front[j].z[1] > front[i].z[1] && front[j].z[2] >= front[i].z[2])) "All solutions in front aren't mutually non-dominated"
        end
    end
end

function DEBUG_two_sols_sorted(sol1::Sol, sol2::Sol)
    @assert sol1.z[1] <= sol2.z[1] "Solutions not sorted wrt z1"
end

function DEBUG_correct_evaluations(population::Vector{Sol}; crowding = false, hv = false)

    pop_copy = copy(population)

    for sol in pop_copy # init values checked
        sol.rank = -1
        if crowding sol.crowding = -1 end
        if hv sol.hv_contrib = -1 end
    end

    F = INDICATOR_fast_non_dominated_sort!(pop_copy)
    for front in F
        crowding && INDICATOR_crowding_distance!(front)
        hv       && INDICATOR_hypervolume_contribution!(front)
    end

    for idx_sol in 1:length(population)
        @assert population[idx_sol].rank == pop_copy[idx_sol].rank "Population not well evaluated (wrong ranks)"
        crowding && @assert population[idx_sol].crowding == pop_copy[idx_sol].crowding "Population not well evaluated (wrong crowding)"
        hv       && @assert population[idx_sol].hv_contrib == pop_copy[idx_sol].hv_contrib "Population not well evaluated (wrong hypervolume contribution)"
    end
end

function DEBUG_unicity(population::Vector{Sol})

    for i in 1:length(population)
        for j in i+1:length(population)
            @assert population[i] != population[j] "No unicity in population"
            @assert population[i].z != population[j].z "No unicity : $(population[i]), $(population[j])"
        end
    end
    
end