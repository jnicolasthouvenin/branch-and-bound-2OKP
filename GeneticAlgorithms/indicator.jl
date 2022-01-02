
#------------- EFFICIENT COMPUTATIONS -----------------#

# The following functions are used in the resolutions,
# thus they do efficient computations.

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

#=

#------------- EXPERIMENTAL COMPUTATIONS -----------------#

# The following functions are used in the experimentations only,
# thus they do expensive and non-optimized computations.

# evaluates the absolute hypervolume of the given front
# `Assertion` : front is a front, hence solutions are mutually non dominated
function INDICATOR_hypervolume(front::Vector{Sol})
    YN = UTIL_front_to_YN(front) # convert into a array of points
    return INDICATOR_hypervolume(YN) # call the hypervolume indicator
end

# evaluates the absolute hypervolume of the given front
# `Assertion` : front is a front, hence solutions are mutually non dominated
function INDICATOR_hypervolume(YN::Vector{Vector{Int}})
    N = length(YN)

    # sort front wrt z1
    perm = sortperm(YN, by = vec -> vec[1])

    hv = YN[perm[1]][1] * YN[perm[1]][2] # initialize the hypervolume

    for idx in 2:N
        hv += (YN[perm[idx]][1]-YN[perm[idx-1]][1])*(YN[perm[idx]][2])
    end

    return hv
end

# indicator denoting the difference between the hypervolume of the reference and the hypervolume of the front
# the smaller the better (like r2 and unary epsilon)
function INDICATOR_hypervolume_inverted(prob::_MOMKP, front::Vector{Sol})
    Y_N = UTIL_read_YN(prob)
    hv_ref = INDICATOR_hypervolume(Y_N)
    hv     = INDICATOR_hypervolume(front)
    return hv_ref  - hv
end

# evaluates the epsilon indicator value of the given front
# `Assertion` : front is a front, hence solutions are mutually non dominated
function INDICATOR_unary_epsilon(prob::_MOMKP, front::Vector{Sol})
    
    debug && DEBUG_feasible_solutions(front)
    debug && DEBUG_is_front(front)

    YN = UTIL_front_to_YN(front) # convert front into a array of points
    ref_YN = UTIL_read_YN(prob) # get reference front

    return INDICATOR_unary_epsilon(YN, ref_YN)

end

# evaluates the epsilon indicator value of the given front
# `Assertion` : YN is a front, hence solutions are mutually non dominated
function INDICATOR_unary_epsilon(YN::Vector{Vector{Int}}, ref_YN::Vector{Vector{Int}})

    N = length(YN)
    N_true = length(ref_YN)

    eps_given_i_j = 0.0 # largest difference given i,j
    eps_given_i = 0.0 # smallest eps_given_i_j, given j
    eps = 0.0 # largest eps_given_i, given i i.e. smallest shift necessary to have YN above ref_YN
    
    buffer = 0.0

    for i in 1:N_true
        #println("-------------------")
        for j in 1:N
            for k in 1:2
                buffer = ref_YN[i][k] - YN[j][k] # maximisation problem

                # update eps_given_i_j (to largest value possible)
                if k == 1 # first objective
                    eps_given_i_j = buffer
                elseif buffer > eps_given_i_j
                    eps_given_i_j = buffer
                end
            end

            # update espJ (to smallest value possible - closest point)
            if j == 1 # first sol of YN
                eps_given_i = eps_given_i_j
            elseif eps_given_i_j < eps_given_i
                eps_given_i = eps_given_i_j
            end
        end

        # update eps (to largest value possible)
        if i == 1 # first sol of ref YN
            eps = eps_given_i
        elseif eps_given_i > eps
            eps = eps_given_i
        end
    end
    
    return eps

end

# evaluates the R2 indicator value of the given front
# `Assertion` : front is a front, hence solutions are mutually non dominated
function INDICATOR_R2(prob::_MOMKP, front::Vector{Sol})

    debug && DEBUG_feasible_solutions(front)
    debug && DEBUG_is_front(front)

    YN = UTIL_front_to_YN(front) # convert front into a array of points
    ref_YN = UTIL_read_YN(prob) # get reference front

    return INDICATOR_R2(YN, ref_YN)

end

# evaluates the R2 indicator value of the given front
# `Assertion` : front is a front, hence solutions are mutually non dominated
function INDICATOR_R2(YN::Vector{Vector{Int}}, ref_YN::Vector{Vector{Int}})

    function generate_lambdas(n::Int)
        lambdas = Matrix{Float64}(undef, n, 2)
        for idx in 1:n
            a = 1.0 * (idx-1) / (n-1)
            lambdas[idx,1] = a
            lambdas[idx,2] = 1 - a
        end
        return lambdas
    end

    L = 10 # number of lambdas
    lambdas = generate_lambdas(L)

    N     = length(YN)
    ref_N = length(ref_YN)

    # invert the front because we want to compute R2 in minimization
    inv_YN = UTIL_invert_front(YN, ref_YN, already_normalized = false)
    inv_ref_YN = UTIL_invert_front(ref_YN, ref_YN, already_normalized = false)

    # for each lambda λ
    # STEP 1 : compute the utility u*(λ,YN) = smallest distance from
    # (0,0) to the any point of inv_YN, using the scalarization vector λ
    # STEP 2 : same of u_ref*(λ, inv_ref_YN)
    # STEP 3 : sum the difference of u - u_ref
    sum = 0.0
    for idx_λ in 1:L

        # compute utility for YN and λ
        u = UTIL_inf()
        for idx_pt in 1:N
            u = min(u, lambdas[idx_λ,1] * inv_YN[idx_pt][1] + lambdas[idx_λ,2] * inv_YN[idx_pt][2])
        end

        # compute utility for ref_YN and λ
        u_ref = UTIL_inf()
        for idx_pt in 1:ref_N
            u_ref = min(u_ref, lambdas[idx_λ,1] * inv_ref_YN[idx_pt][1] + lambdas[idx_λ,2] * inv_ref_YN[idx_pt][2])
        end

        @assert u >= u_ref "Impossible ! $u > $u_ref"

        sum += u - u_ref

    end

    return sum / L # return the normalized value of the sum

end

# evalutes the front with regards to the three considered quality indicators (hypervolume, unary epsilon, R2)
function INDICATOR_evaluate_front_estimation(prob::_MOMKP, front::Vector{Sol})

    debug && DEBUG_feasible_solutions(front)
    debug && DEBUG_is_front(front)

    #hv  = INDICATOR_hypervolume_inverted(prob, front)
    hv  = 1 - (INDICATOR_hypervolume(front) / INDICATOR_hypervolume(UTIL_read_YN(prob)))
    eps = INDICATOR_unary_epsilon(prob, front)
    r2  = INDICATOR_R2(prob, front)

    return hv, eps, r2

end

=#