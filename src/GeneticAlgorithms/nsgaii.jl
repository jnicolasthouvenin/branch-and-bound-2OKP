
function NSGAII_update(P::Vector{Sol}, config::Config)

    # P is supposed to be already evaluated (vector of AVLs)
    debug && DEBUG_correct_evaluations(P, crowding = true)

    #------ Offsprings -------#

    N = length(P)

    # generate offsprings
    offsprings = Vector{Sol}(undef, N)
    idx_offspring = 0
    while idx_offspring < N

        # choose parents
        parent_1 = SELECTION_binary_tournament(P)
        parent_2 = SELECTION_binary_tournament(P)
        while parent_2 == parent_1
            parent_2 = SELECTION_binary_tournament(P)
        end

        # recombine parents
        offspring = config.crossover(parent_1, parent_2)

        # mutate offspring
        if rand() < config.mutation_rate
            MUTATION_flip!(offspring)
        end

        # repair offspring
        config.reparation!(offspring)

        # check that the offspring is not already in the population
        is_new = true
        for sol in P
            if sol.z == offspring.z
                is_new = false
                break
            end
        end

        if is_new
            # store offspring
            idx_offspring += 1
            offsprings[idx_offspring] = offspring
        end
    end

    #------ Create merged population R, sort it and evaluate it ------#

    R = UTIL_concatenate_solutions_vectors(P,offsprings)

    F = INDICATOR_fast_non_dominated_sort!(R)

    min_size_new_pop = min(2*N, Int(floor(length(F[1])*1.5))) # compute new adaptative size

    N = max(N, min_size_new_pop)

    new_P = Vector{Sol}(undef, N)
    size_new_P = 0
    idx_front = 1

    while size_new_P + length(F[idx_front]) <= N # we can the whole front

        # evaluate front (for next iteration)
        INDICATOR_crowding_distance!(F[idx_front])

        # add the whole front
        for sol_front in F[idx_front]
            size_new_P += 1
            new_P[size_new_P] = sol_front
        end

        idx_front += 1
    end

    if size_new_P < N # we still have space for a few solutions of the next front
        # evaluate front
        INDICATOR_crowding_distance!(F[idx_front])

        # sort front by decreasing crowding distance
        perm = sortperm(F[idx_front], by = sol -> sol.crowding, rev = true)

        idx_sol = 0
        while size_new_P < N
            size_new_P += 1
            idx_sol += 1
            new_P[size_new_P] = F[idx_front][perm[idx_sol]]
        end

    end

    debug && DEBUG_correct_evaluations(new_P, crowding = true)

    return new_P

end

function NSGAII_solve(prob::_MOMKP, config::Config; seeding_sols = Vector{Sol}(undef, 0), max_t = 10.0, verbose = true)

    st = time()

    # create population
    P = INIT_random_feasible_population(config.size_pop_init, prob, seeding_sols = seeding_sols)

    # evaluate population
    F = INDICATOR_fast_non_dominated_sort!(P)
    for front in F
        INDICATOR_crowding_distance!(front)
    end

    # iterate
    it = 0
    while !UTIL_termination_criteria_met(time()-st, max_t)
        it += 1
        P = NSGAII_update(P, config)
    end

    # evaluate final front
    F = INDICATOR_fast_non_dominated_sort!(P)
    INDICATOR_crowding_distance!(F[1])

    last_front = UTIL_remove_doublons_from_front(F[1])
    sort!(last_front, by = sol -> sol.z[1])

    verbose && println("total time spent : ",time() - st)
    verbose && println("total nb of iter : ",it)
    verbose && println("total nb of sols : ",length(last_front),"\n")

    return last_front

end