
""" Functions that write and read 2OKP instances on files """

# All functions are public

"""
Returns the name of the file associated with the given problem `prob` and `id`

`id`, unique number to avoid writing on existing files (generally between 1 and 5)
"""
function FILES_getProbName(prob::BiOKP,id::Int;isInteger=true)
    if isInteger
        return string("Inst_",string(prob.nbVar),"_",string(id))
    else
        return string("Inst_",string(prob.nbVar),"_",string(id),"_float")
    end
end

"""
Write the given problem `prob` attributes on a file

`id`, unique number to avoid writing on existing files (generally between 1 and 5)

# Example
```jldoctest
julia> FILES_writeProblem(myProblem,10)
```
"""
function FILES_writeProblem(prob::BiOKP,id::Int; isInteger = true)
    fileName = FILES_getProbName(prob,id,isInteger=isInteger)

    f = open(string("instances/",fileName,".dat"),"w")

    println(f,string(prob.nbVar," ",prob.nbObj))
    for i in 1:2
        for profit in prob.profits[i,:]
            print(f,string(profit," "))
        end
        println(f,"")
    end
    for weight in prob.weights
        print(f,string(weight," "))
    end
    println(f,"")
    println(f,string(prob.maxWeight))
    println(f,string(prob.isInteger))

    close(f)
end

"""
Returns the problem stored in the given file

The `fileName` is the name of the given file (no extension, no path, just the name)

# Example
```jldoctest
julia> p = FILES_readProblem("Inst_5_1")
BiOKP(5, 2, [29.0 28.0 … 18.0 31.0; 11.0 36.0 … 27.0 37.0], [16.0, 50.0, 32.0, 13.0, 50.0], 65.0, true)
```
"""
function FILES_readProblem(fileName::String)

    f = open(string("instances/",fileName,".dat"),"r")

    first_line::Array{SubString{String},1} = split(readline(f))
    nbVars = parse(Int64, first_line[1])
    nbObj = parse(Int64, first_line[2])
    profits_1 = split(readline(f))
    profits_2 = split(readline(f))
    profits = zeros(Float64, 2, nbVars)
    weights_str = split(readline(f))
    weights = zeros(Float64,nbVars)
    for i in 1:nbVars
        profits[1,i] = parse(Float64,profits_1[i])
        profits[2,i] = parse(Float64,profits_2[i])
        weights[i] = parse(Float64,weights_str[i])
    end
    maxWeight = parse(Float64,readline(f))
    isInteger = parse(Bool,readline(f))

    return BiOKP(nbVars,nbObj,profits,weights,maxWeight,isInteger)
end