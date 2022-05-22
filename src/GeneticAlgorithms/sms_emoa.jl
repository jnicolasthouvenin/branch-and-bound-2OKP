
# SMS EMOA stands for S Metric Selection Evolutionary Multi-Objective Algorithm
# S metric is the measure of hypervolume
# Simply put, each iteration, the algorithm generate a random feasible offspring and insert it in the popultion,
# and then removes the individual that is the less contributing to the overall hypervolume.
# Thus: |P_k| = |P_{k+1}| (stationary population)
# The algorithm tend to focus more on the "eye" of the pareto front.
# source: https://doi.org/10.1007/978-3-540-31880-4_5

function SMS_EMOA_update(P::Vector{Sol}, config::Config)

    debug && DEBUG_feasible_solutions(P)
    debug && DEBUG_correct_evaluations(P, hv = true)

    N = length(P)

    #--------- generate offspring ----------#

    offspring_has_been_generated = false
    offspring = 0
    while !offspring_has_been_generated

        # choose parents
        parent_1 = SELECTION_binary_tournament(P, criteria = HV)
        parent_2 = SELECTION_binary_tournament(P, criteria = HV)
        while parent_2 == parent_1
            parent_2 = SELECTION_binary_tournament(P, criteria = HV)
        end

        # recombine parents
        offspring = config.crossover(parent_1, parent_2)

        # mutate offspring
        if rand() < config.mutation_rate
            MUTATION_flip!(offspring)
        end

        # repair offspring
        config.reparation!(offspring)

        # check if the offspring is new to the population
        offspring_has_been_generated = true
        for sol in P
            if sol.z == offspring.z
                offspring_has_been_generated = false
                break
            end
        end
    end

    #-------------- add offspring to current population ---------------#

    # add offspring to the population
    push!(P, offspring)

    #- sort population by rank, and compute hypervolume for each rank -#

    # evaluate population
    F = INDICATOR_fast_non_dominated_sort!(P)
    for front in F
        INDICATOR_hypervolume_contribution!(front)
    end

    R_I = F[length(F)]

    # find the solution in R_I that contributes the less to the hypervolume
    idx_min = 1
    sol_min = R_I[1]
    hv_min = R_I[1].hv_contrib
    for idx in 1:length(R_I)
        sol = R_I[idx]
        if sol.hv_contrib < hv_min
            idx_min = idx
            sol_min = sol
            hv_min = sol.hv_contrib
        end
    end

    deleteat!(R_I, idx_min) # remove the worst contributing solution in the front

    # evaluate the new contributions (for next iteration)
    if length(R_I) > 0
        INDICATOR_hypervolume_contribution!(R_I)
    end

    # mark the solution that contributes the less to the hypervolume
    sol_min.hv_contrib = -1

    # find the solution marked in P
    idx_worst_contribution_sol = 0
    idx = 1
    while idx <= N+1 && idx_worst_contribution_sol == 0
        if P[idx].hv_contrib == -1 # we found the marked solution
            idx_worst_contribution_sol = idx
        end
        idx += 1
    end

    @assert idx_worst_contribution_sol != 0 "The marked solution isn't found"

    #------------- Remove the worst contribution solution from P -----------#

    # remove the marked solution (contributing less to the hypervolume)
    deleteat!(P, idx_worst_contribution_sol)

    debug && DEBUG_feasible_solutions(P)
    debug && DEBUG_correct_evaluations(P, hv = true)

    return P
end

function SMS_EMOA_solve(prob::_MOMKP, config::Config; seeding_sols = Vector{Sol}(undef, 0), max_t = 10.0, verbose = true)

    st = time()

    # create population
    P = INIT_random_feasible_population(config.size_pop_init, prob, seeding_sols = seeding_sols)

    # evaluate population
    F = INDICATOR_fast_non_dominated_sort!(P)
    for front in F
        INDICATOR_hypervolume_contribution!(front)
    end

    # iterate
    it = 0
    while !UTIL_termination_criteria_met(time()-st, max_t)
        it += 1
        P = SMS_EMOA_update(P, config)
    end

    # evaluate population
    F = INDICATOR_fast_non_dominated_sort!(P)

    last_front = UTIL_remove_doublons_from_front(F[1])

    sort!(last_front, by = sol -> sol.z[1])

    verbose && println("total time spent : ",time() - st)
    verbose && println("total nb of iter : ",it)
    verbose && println("total nb of sols : ",length(last_front),"\n")

    return last_front
end