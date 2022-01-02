
@enum Criteria CROWDING HV # different criterium for SELECTION_binary_tournament

struct _MOMKP
    name :: String # name of the instance
    n  :: Int # number of items
    m  :: Int # number of dimensions
    P  :: Matrix{Int} # profit of items for the objectives, k=1..p, j=1..n
    W  :: Matrix{Int} # weight of items for the constraints, i=1..m, j=1..n
    ω  :: Vector{Int} # capacity of knapsacks, i=1..m
end

#=
struct _bi01IP
    C  :: Matrix{Int} # objective functions, k=1..2, j=1..n
    A  :: Matrix{Int} # matrix of constraints, i=1..m, j=1..n
    b  :: Vector{Int} # right-hand side, i=1..m
end
=#

#----------- Constructors --------------#

#=
function MOMKP(instance_file_path::String)

    if instance_file_path[1:11] == "instancesPG"
        momkp = readInstanceMOMKPformatPG(instance_file_path)
    elseif instance_file_path[1:11] == "instancesZL"
        momkp = readInstanceMOMKPformatZL(instance_file_path)
    else
        @assert false "Incorrect instance file path : $instance_file_path"
    end

    return momkp
end=#

#=
function bi01KP(instance_file_path::String)
    if instance_file_path[1:11] == "instancesPG"
        momkp = readInstanceMOMKPformatPG(instance_file_path)
        bi01kp = _bi01IP(momkp.P, momkp.W, momkp.ω)
    elseif instance_file_path[1:11] == "instancesZL"
        momkp = readInstanceMOMKPformatZL(instance_file_path)
        bi01kp = _bi01IP(momkp.P, momkp.W, momkp.ω)
    else
        @assert false "Incorrect instance file path : $instance_file_path"
    end

    return bi01kp
end=#

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

#--------------- Parsers -----------------#

#=
# Credit : Xavier Gandibleux - https://github.com/xgandibleux
function readInstanceMOMKPformatZL(fname::String)

    f = open(fname)
    lines = readlines(f)
    close(f)

    # Extract the problem dimensions from the file name ------------------------
    line=(split(fname,"."))
    #line = split(lines[1]," ")
    m = parse(Int, line[3])
    p = m
    n = parse(Int, line[2])
    momkp = _MOMKP(fname,n,m,zeros(Int,p,n),zeros(Int,m,n),zeros(Int,m))

    # Extract the data of the MKP from the following lines ---------------------
    i=0; j=0
    for l = 2:length(lines)
        if occursin("knapsack ", lines[l])
            i = i + 1
            j = 0
        elseif occursin("capacity", lines[l])
            line = split(lines[l],"+")
            momkp.ω[i] = floor(Int, parse(Float64, line[2]))
        elseif occursin("item", lines[l])
            j = j + 1
        elseif occursin("weight", lines[l])
            line = split(lines[l],"+")
            momkp.W[i,j] = parse(Int, line[2])
        elseif occursin("profit", lines[l])
            line = split(lines[l],"+")
            momkp.P[i,j] = parse(Int, line[2])
        end
    end

    return momkp
end

# Credit : Xavier Gandibleux - https://github.com/xgandibleux
function readInstanceMOMKPformatPG(fname::String)

    f=open(fname)

    n = parse(Int, readline(f) )
    m = 2  # 2 constraints
    p = 2  # 2 objectives

    momkp = _MOMKP(fname,n,m,zeros(Int,p,n),zeros(Int,m,n),zeros(Int,m))

    momkp.ω[1] = parse(Int, readline(f) )
    momkp.ω[2] = parse(Int, readline(f) )

    # Extract the data of the MKP from the following lines ---------------------
    for j = 1:n
        momkp.P[1,j],
        momkp.P[2,j],
        momkp.W[1,j],
        momkp.W[2,j] = parse.(Int, split(readline(f)) )
    end
    close(f)

    return momkp
end=#