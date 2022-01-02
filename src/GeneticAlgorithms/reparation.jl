
function REPARATION_basic!(offspring::Sol)

    debug && DEBUG_correct_solution(offspring)

    if !SOLUTION_is_feasible(offspring)
        # repair solution
        for var in 1:offspring.prob.n
            if offspring.x[var] == 0
                if SOLUTION_is_feasible_after_adding_item(offspring, var)
                    SOLUTION_add_item!(offspring, var)
                end
            else
                if !SOLUTION_is_feasible(offspring)
                    SOLUTION_remove_item!(offspring, var)
                end
            end
        end
    end

    debug && DEBUG_feasible_solution(offspring)
    
end

function REPARATION_advanced!(offspring::Sol)

    debug && DEBUG_correct_solution(offspring)

    if !SOLUTION_is_feasible(offspring)
        
        λ = rand() # we pick a random direction for the reparation
        momkp = offspring.prob

        # count the number of ones (thus also the number of zeros) in the solution
        nb_ones = 0
        for component in offspring.x
            component == 1 ? nb_ones += 1 : nothing
        end

        idx_of_ones = Vector{Int}(undef, nb_ones) # vars to one
        idx_ones = 0

        idx_of_zeros = Vector{Int}(undef, momkp.n-nb_ones) # vars to zero
        idx_zeros = 0

        # fill array idx_of_ones and idx_of_zeros
        for var in 1:momkp.n
            if offspring.x[var] == 1
                idx_ones += 1
                idx_of_ones[idx_ones] = var
            else
                idx_zeros += 1
                idx_of_zeros[idx_zeros] = var
            end
        end

        # compute the ratio of each variable : ratios[idx_var]
        ratios = Vector{Float64}(undef, momkp.n)
        for var in 1:momkp.n
            profit = λ*momkp.P[1, var] + (1-λ)*momkp.P[2, var]
            weight = 0
            for dim in 1:momkp.m
                weight += momkp.W[dim, var]
            end
            ratios[var] = profit/weight
        end

        # sort the ratios of ones (worst to best)
        sort!(idx_of_ones, by = idx_var -> ratios[idx_var])

        # sort the ratios of zeros (best to worst)
        sort!(idx_of_zeros, by = idx_var -> ratios[idx_var], rev = true)

        # remove items to repair the solution starting from the worst variables
        idx = 0
        while idx < length(idx_of_ones) && !SOLUTION_is_feasible(offspring)
            idx += 1
            SOLUTION_remove_item!(offspring, idx_of_ones[idx])
        end

        # add items in the solution (while remaining feasible) starting from the best variables
        idx = 1
        while idx <= length(idx_of_zeros) && SOLUTION_is_feasible_after_adding_item(offspring, idx_of_zeros[idx])
            SOLUTION_add_item!(offspring, idx_of_zeros[idx])
            idx += 1
        end
    end

    debug && DEBUG_feasible_solution(offspring)
    
end