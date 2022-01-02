
""" Declaration of the original types used in the algorithm """

################################ ENUMS ################################

@enum Pruned OPTIMALITY INFEASIBILITY DOMINANCE NONE

################################ STRUCTS ################################

# for storing a problem
struct BiOKP
	nbVar::Int
	nbObj::Int
	profits::Array{Float64,2}
	weights::Array{Float64,1}
	maxWeight::Float64
	isInteger::Bool
end

# for storing a solution
mutable struct Sol
	x::Vector{Float64}
	y::Vector{Float64}
	w::Float64
	isBinary::Bool
end

# for storing a specific assignment for a BiOKP
mutable struct Assignment
	assign::Vector{Int}
	profit::Vector{Float64}
	weight::Float64
	lastAssigned::Int
end

# for storing a pair of solutions
struct PairOfSolution
	solL::Sol
	solR::Sol
end

# for storing upper bounds
struct DualSet
	A::Array{Float64, 2}
	b::Array{Float64, 1}
end

# used in the parametric pretreatment
mutable struct LambdaChange
	index::Tuple{Int, Int}
	value::Frac
end

# for storing parametric pretreatment
mutable struct LinearSpeedUp
    indicies::Vector{Vector{Tuple{Int,Int}}}
    permObj::Vector{Int}
end

# for storing known solutions that are present in the parent and the child node
mutable struct ParentToChild
    knownSols::LinkedList{Sol} # solutions that are identical from the parent to the child node
    previousIndicies::LinkedList{Int} # indicies of the knownSols inside the parent UB
    rightSol::Sol # the last solution of knownSols list. The one on the right.
    leftKnown::Bool # the left sol of the parent is a known sol
    rightKnown::Bool # the right sol of the parent is a known sol
end

# for couting the number of subproblems studied
mutable struct Iterator
	value::Int
end

################################ LOWER BOUND  ################################

mutable struct LowerBound
    sols::LinkedList{Sol}
end

################################ UPPER BOUND  ################################

mutable struct UpperBound
    sols::LinkedList{Sol}
    segments::DualSet
    explicit::Bool # wether we use the explicit fathoming test of not, if false, then segments is not defined
end

################################ CONSTRUCTORS ################################

# PROBLEM ######################

function BiOKP()
    return BiOKP(6, 2, [11 2 8 10 9 1; 2 7 8 4 1 3], [4, 4, 6, 4, 3, 2], 11, true)
end

function BiOKP(nbVar::Int, nbObj::Int, profits::Array{Float64,2}, weights::Array{Float64, 1}, maxWeight::Float64)
	return BiOKP(nbVar, nbObj, profits, weights, maxWeight, false)
end

function BiOKP(fname::String)
	f = open(fname)

    nbVar, nbObj = parse.(Int, split(readline(f)))

	profits = zeros(Float64, nbObj, nbVar)
	weights = zeros(Float64, nbVar)

	for iter = 1:nbObj
		profits[iter, 1:end] = parse.(Float64, split(readline(f)))
	end

	weights = parse.(Float64, split(readline(f)))

	maxWeight = parse(Float64, readline(f))

    readline(f) == "true" ? isInteger = true : isInteger = false

	close(f)

	return BiOKP(nbVar, nbObj, profits, weights, maxWeight, isInteger)
end

function BiOKP(nbVar::Int; isInteger::Bool = true)
    if isInteger
        weights = rand(1:50, nbVar)
        maxWeight = Int(ceil(sum(weights) / 2.5))
        return BiOKP(nbVar, 2, rand(1:50, 2, nbVar), weights, maxWeight <= 50 ? 51 : maxWeight, true)
    else
        weights = rand(Float64, nbVar) * 49 .+ 1
        maxWeight = sum(weights) / 2.5
        return BiOKP(nbVar, 2, rand(Float64, 2, nbVar) * 49 .+ 1, weights, maxWeight <= 50 ? 51 : maxWeight, false)
    end
end

# SOLUTION ######################

function Sol()
	return Sol( Vector{Float64}(), Vector{Float64}(), 0., false)
end

function Sol(x::Vector{Float64}, y::Vector{Float64}, w::Float64)
	return Sol(x, y, w, true)
end

# ASSIGNMENT ######################

function Assignment()
	return Assignment(Vector{Int}(), Vector{Float64}(), 0, 0)
end

function Assignment(prob::BiOKP)
	return Assignment(ones(Int, prob.nbVar)*-1, zeros(Float64, prob.nbObj), 0, 0)
end

# PAIR OF SOLUTION ######################

function PairOfSolution()
	return PairOfSolution(Sol(),Sol())
end

# DUAL SET #############################

function DualSet()
    return DualSet(Matrix{Float64}(undef, 0, 2), Vector{Float64}(undef, 0))
end

# LINEAR SPEED UP ######################

function LinearSpeedUp()
    return LinearSpeedUp(Vector{LambdaChange}(),Vector{Int}())
end

# FROM PARENT TO CHILD ######################

function ParentToChild()
    return ParentToChild(nil(Sol),nil(Int),Sol(),false,false)
end

# ITERATOR ######################

function Iterator()
	return Iterator(0)
end

# LOWER BOUND ######################

function LowerBound()
    return LowerBound(nil(Sol))
end

# UPPER BOUND ######################

function UpperBound()
    return UpperBound(nil(Sol), DualSet(), false)
end

function UpperBound(sols::LinkedList{Sol})
    return UpperBound(sols, DualSet(), false)
end

function UB_setSegments!(UB::UpperBound, segments::DualSet)
    UB.segments = segments
    UB.explicit = true
end

################################ OPERATORS ################################

import Base.:(==), Base.:(!=)

Base.:(==)(sol1::Sol, sol2::Sol) = sol1.y == sol2.y
Base.:(==)(pair1::PairOfSolution, pair2::PairOfSolution) = pair1.solL == pair2.solL && pair1.solR == pair2.solR
Base.:(==)(x::DualSet, y::DualSet) = (x.A == y.A) && (x.b == y.b)
Base.:(!=)(pair1::PairOfSolution, pair2::PairOfSolution) = !(pair1 == pair2)

################################ COPY ################################

import Base.copy
Base.copy(sol::Sol) = Sol(copy(sol.x), copy(sol.y), sol.w, sol.isBinary)

################################ SHOW ################################

# ELEMENTS ######################

function Base.show(io::IO, sol::Sol)
    index = 0
    i = 1
    while index == 0 && i < length(sol.x)
        if sol.x[i] > 0. && sol.x[i] < 1.
            index = i
        end
        i += 1
    end
    print(io,"$(sol.y)")
end

function Base.show(io::IO, pair::PairOfSolution)
    print(io,"($(pair.solL.y),$(pair.solR.y))")
end

# LISTS ######################

function Base.show(io::IO, linkedL::LinkedList{Sol})
    print(io, "\n")
    for sol in linkedL
        println(io, "$(sol)")
    end
    println(io, "")
end

function Base.show(io::IO, vec::Vector{Sol})
    print(io, "\n")
    for sol in vec
        println(io, "$(sol)")
    end
    println(io, "")
end

function Base.show(io::IO, vec::Vector{PairOfSolution})
    print(io,"\n")
    for elt in vec
        println(io, "$(elt.solL.y),$(elt.solR.y)")
    end
    println(io,"")
end

function Base.show(io::IO, vec::Vector{LambdaChange})
    print(io,"(")
    for elt in vec
        println(io, "$(elt.value) - $(elt.index)")
    end
    println(io,")")
end