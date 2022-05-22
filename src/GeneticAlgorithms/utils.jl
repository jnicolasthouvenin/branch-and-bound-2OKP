
UTIL_inf() = 2^31 - 1

function UTIL_concatenate_solutions_vectors(population::Vector{Sol}, offsprings::Vector{Sol})

    debug && DEBUG_feasible_solutions(population)
    debug && DEBUG_feasible_solutions(offsprings)

    result = vcat(population,offsprings)

    return result
end

#=
function UTIL_normalize_front(YN::Vector{Vector{Int}}, ref_YN::Vector{Vector{Int}})

    N = length(YN)
    normalized_YN = Vector{Vector{Float64}}(undef, N)

    #minimums = [ref_YN[1][1]              - 1 , ref_YN[length(ref_YN)][2] - 1 ]
    maximums = [ref_YN[length(ref_YN)][1] + 1 , ref_YN[1][2]              + 1 ]

    for idx in 1:N
        normalized_YN[idx] = [0.0, 0.0] # init
        normalized_YN[idx][1] = YN[idx][1] / maximums[1]
        normalized_YN[idx][2] = YN[idx][2] / maximums[2]
    end

    #println("normalized_YN = ",normalized_YN)

    return normalized_YN
end=#

#=
function UTIL_invert_front(YN::Vector{Vector{Int}}, ref_YN::Vector{Vector{Int}}; already_normalized = true)

    if !already_normalized
        YN = UTIL_normalize_front(YN, ref_YN)
    end

    N = length(YN)
    inversed_YN = YN[1:N] # deepcopy

    for idx in 1:N
        @assert YN[idx][1] >= 0 && YN[idx][1] <= 1 "Front is not normalized"
        @assert YN[idx][2] >= 0 && YN[idx][2] <= 1 "Front is not normalized"

        inversed_YN[idx][1] = 1 - inversed_YN[idx][1]
        inversed_YN[idx][2] = 1 - inversed_YN[idx][2]
    end

    sort!(inversed_YN, by = vec -> vec[1])

    return inversed_YN
end
=#

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

#=
function UTIL_front_to_YN(front::Vector{Sol})

    # convert front into a array of points
    YN = Vector{Vector{Int}}(undef, length(front))
    for i in 1:length(front)
        YN[i] = [front[i].z[1], front[i].z[2]]
    end

    return YN

end

function UTIL_write_YN(prob::_MOMKP, YN::Vector{Vector{Int}})

    f = open(string("YN/",prob.name),"w")

    println("   write in ",string("YN/",prob.name))

    println(f,string(length(YN)))
    for i in 1:length(YN)
        println(f,string(YN[i][1]," ",YN[i][2]))
    end

    close(f)
end

function UTIL_read_YN(prob::_MOMKP)

    f = open(string("YN/",prob.name),"r")
    lines = readlines(f)
    close(f)

    len = length(lines)-1

    YN = Vector{Vector{Int}}(undef, len)

    for i in 1:len
        YN[i] = parse.(Int64, split(lines[i+1]))
    end

    return YN

end
=#