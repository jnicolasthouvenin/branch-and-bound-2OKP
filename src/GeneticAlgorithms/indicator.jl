
# Functions for computing indicators on solutions

# true is sol1 dominates sol2
function INDICATOR_is_dominating(sol1::Sol, sol2::Sol)

    debug && DEBUG_correct_solutions([sol1,sol2])

    return (sol1.z[1] >= sol2.z[1] && sol1.z[2] > sol2.z[2]) || (sol1.z[1] > sol2.z[1] && sol1.z[2] >= sol2.z[2])
end

# evaluates the rank of each solution in the given population
function INDICATOR_fast_non_dominated_sort!(population::Vector{Sol})

    debug && DEBUG_feasible_solutions(population)

    # true if no individual in population has rank idx_rank
    function rank_empty(population::Vector{Sol}, idx_rank::Int)
        for p in population
            if p.rank == idx_rank
                return false
            end
        end
        return true
    end

    # initialize the rank of each sol in the population
    for p in population
        p.rank = 0
    end

    n = length(population)
    S = [ [] for idx in 1:n] # S[i] list of indices of sols dominated by pop[i]
    domination_count = zeros(Int, n) # ...[i] number of sols that dominates pop[i]

    for idx_p in 1:n
        p = population[idx_p]
        for idx_q in idx_p+1:n
            q = population[idx_q]
            if INDICATOR_is_dominating(p, q) || p.z == q.z
                push!(S[idx_p], idx_q)
                domination_count[idx_q] += 1
            elseif INDICATOR_is_dominating(q, p)
                push!(S[idx_q], idx_p)
                domination_count[idx_p] += 1
            end
        end

        if domination_count[idx_p] == 0
            p.rank = 1
        end
    end

    idx_rank = 1

    while !rank_empty(population, idx_rank)
        for idx_p in 1:n
            if population[idx_p].rank == idx_rank
                for idx_q in S[idx_p]
                    domination_count[idx_q] -= 1
                    if domination_count[idx_q] == 0
                        population[idx_q].rank = idx_rank + 1
                    end
                end
            end
        end
        idx_rank += 1
    end

    # compute explicite fronts
    nb_ranks = idx_rank - 1
    F = Vector{Vector{Sol}}(undef, nb_ranks)
    for idx in 1:nb_ranks
        F[idx] = Vector{Sol}(undef, 0)
        for p in population
            if p.rank == idx
                push!(F[idx], p)
            end
        end
    end

    return F

end

# evaluates the crowding distance of each solution in the given front
# `Assertion` : front is a front, hence solutions are mutually non dominated
function INDICATOR_crowding_distance!(front::Vector{Sol})

    debug && DEBUG_feasible_solutions(front)
    debug && DEBUG_is_front(front)

    n = length(front)

    if n == 1

        front[1].crowding = UTIL_inf()

    else

        # initialize crowding distance for each solution
        for sol in front
            sol.crowding = 0
        end

        # permutation to sort front with regards to the first objective
        perm = sortperm(front, by = sol -> sol.z[1]) # O(nlog(n)))

        # set crowding of extreme sols to Infinity to always select them
        front[perm[1]].crowding = UTIL_inf()
        front[perm[n]].crowding = UTIL_inf()

        # first objective
        for idx in 2:n-1
            front[perm[idx]].crowding += front[perm[idx+1]].z[1] - front[perm[idx-1]].z[1]
        end

        # second objective
        for idx in 2:n-1
            front[perm[n+1-idx]].crowding += front[perm[n+1-(idx+1)]].z[2] - front[perm[n+1-(idx-1)]].z[2]
        end

    end

end

# evaluates the hypervolume contribution of each solution in the given front
# `Assertion` : front is a front, hence solutions are mutually non dominated
function INDICATOR_hypervolume_contribution!(front::Vector{Sol})

    debug && DEBUG_feasible_solutions(front)
    debug && DEBUG_is_front(front)

    N = length(front)

    # sort front wrt z1
    perm = sortperm(front, by = sol -> sol.z[1])

    # extreme solutions get an infinite hypervolume contribution
    front[perm[1]].hv_contrib = UTIL_inf()
    front[perm[N]].hv_contrib = UTIL_inf()

    for i in 2:N-1
        front[perm[i]].hv_contrib = (front[perm[i]].z[1] - front[perm[i-1]].z[1]) + (front[perm[i]].z[2] - front[perm[i+1]].z[2])
    end

end
