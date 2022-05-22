
UTIL_inf() = 2^31 - 1

function UTIL_concatenate_solutions_vectors(population::Vector{Sol}, offsprings::Vector{Sol})

    debug && DEBUG_feasible_solutions(population)
    debug && DEBUG_feasible_solutions(offsprings)

    result = vcat(population,offsprings)

    return result
end

function UTIL_remove_doublons_from_front(front::Vector{Sol})

    debug && DEBUG_is_front(front)

    new_front = Vector{Sol}(undef, 0)

    for sol in front
        if !(sol in new_front)
            push!(new_front, sol)
        else
            println("REMOVE DOUBLON")
        end
    end

    debug && DEBUG_is_front(new_front)
    debug && DEBUG_unicity(new_front)

    return new_front
end

function UTIL_termination_criteria_met(it::Int, nb_it::Int)
    return it > nb_it
end

function UTIL_termination_criteria_met(t::Float64, max_t::Float64)
    return t > max_t
end
