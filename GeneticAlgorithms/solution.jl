
function SOLUTION_is_feasible(sol::Sol)

    debug && DEBUG_correct_solution(sol)

    feasible = true
    dim = 1
    while feasible && dim <= sol.prob.m
        if sol.w[dim] > sol.prob.ω[dim]
            feasible = false
        end
        dim += 1
    end
    
    return feasible
end

function SOLUTION_is_feasible_after_adding_item(sol::Sol, idx_var::Int)

    debug && DEBUG_correct_solution(sol)

    @assert sol.x[idx_var] == 0 "The item is already in the bag"
    
    feasible = true
    dim = 1
    while feasible && dim <= sol.prob.m
        if sol.w[dim] + sol.prob.W[dim, idx_var] > sol.prob.ω[dim]
            feasible = false
        end
        dim += 1
    end
    
    return feasible
end

function SOLUTION_add_item!(sol::Sol, item::Int)

    debug && @assert sol.x[item] == 0 "Item already in the bag"
    debug && DEBUG_correct_solution(sol) # abusive

    sol.x[item] = 1
    sol.z[1] += sol.prob.P[1,item]
    sol.z[2] += sol.prob.P[2,item]
    for dim in 1:sol.prob.m
        sol.w[dim] += sol.prob.W[dim,item]
    end

    debug && DEBUG_correct_solution(sol)

end

function SOLUTION_remove_item!(sol::Sol, item::Int)

    debug && @assert sol.x[item] == 1 "Item not in the bag"
    debug && DEBUG_correct_solution(sol) # abusive

    momkp = sol.prob

    sol.x[item] = 0
    sol.z[1] -= momkp.P[1,item]
    sol.z[2] -= momkp.P[2,item]
    for dim in 1:momkp.m
        sol.w[dim] -= momkp.W[dim,item]
    end

    debug && DEBUG_correct_solution(sol)
    
end

function SOLUTION_empty_solution(momkp::_MOMKP)

    empty_sol = Sol(momkp, zeros(Int, momkp.n), [0, 0], zeros(Int, momkp.m))

    return empty_sol
end

function SOLUTION_list_of_empty_solutions(nb_solutions::Int, momkp::_MOMKP)

    vec = Vector{Sol}(undef, nb_solutions)
    for i in 1:nb_solutions
        vec[i] = SOLUTION_empty_solution(momkp)
    end

    return vec
end