
@enum Criteria CROWDING HV # different criterium for SELECTION_binary_tournament

struct _MOMKP
    name :: String # name of the instance
    n  :: Int # number of items
    m  :: Int # number of dimensions
    P  :: Matrix{Int} # profit of items for the objectives, k=1..p, j=1..n
    W  :: Matrix{Int} # weight of items for the constraints, i=1..m, j=1..n
    ω  :: Vector{Int} # capacity of knapsacks, i=1..m
end

mutable struct Sol
    prob::_MOMKP

    x::Vector{Int} # length(x) = prob.n
    z::Vector{Int} # length(z) = 2
    w::Vector{Int}

    rank::Int
    crowding::Int
    hv_contrib::Int

    Sol(prob, x, z, w) = new(prob, x, z, w, 0, 0, 0) # only aload constructor
end

import Base.:(==), Base.:(!=)

Base.:(==)(s1::Sol, s2::Sol) = (s1.x == s1.x) && (s1.z == s2.z)
Base.:(!=)(s1::Sol, s2::Sol) = (s1.x != s1.x) || (s1.z != s2.z)

#------------ Configurations -------------#

struct Config
    size_pop_init::Int
    crossover::Function
    mutation_rate::Float64
    reparation!::Function
end
