
"""
return a parent winner of a binary tournament selection

`Assertion` : the given population must have been evaluated for crowding distance
"""
function SELECTION_binary_tournament(population::Vector{Sol}; criteria = CROWDING)

    debug && DEBUG_feasible_solutions(population)

    function winner_match(population::Vector{Sol}, ind_1::Int, ind_2::Int)
        if INDICATOR_is_dominating(population[ind_1], population[ind_2])
            return ind_1
        elseif INDICATOR_is_dominating(population[ind_2], population[ind_1])
            return ind_2
        elseif criteria == CROWDING
            if population[ind_1].crowding > population[ind_2].crowding
                return ind_1
            else
                return ind_2
            end
        elseif criteria == HV
            if population[ind_1].hv_contrib > population[ind_2].hv_contrib
                return ind_1
            else
                return ind_2
            end
        end
    end

    N = length(population)

    # pick two random distinct solutions
    ind_1 = rand(1:N)
    ind_2 = rand(1:N)
    while ind_1 == ind_2
        ind_2 = rand(1:N)
    end

    # return only the winner
    selected_solution = population[winner_match(population, ind_1, ind_2)]

    return selected_solution

end
