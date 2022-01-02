
""" Implementation of the parametric method for linear relaxation """

# Public functions :
# - initialiseLinearSpeedUp : initialisation of the parametric method
# - parametricMethod : use a LinearSpeedUp object and Assignment to compute a relaxed UB

#-------------------------------------------------------------------#
#------------------------ PUBLIC FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

" Returns a `LinearSpeedUp` object initialised "
function initialiseLinearSpeedUp(prob::BiOKP)

    @timeit to "initialiseLinearSpeedUp" begin

    # compute the ratios
    ratios = Matrix{Frac}(undef,2,prob.nbVar)
    for indexObj in 1:2
        for indexVar in 1:prob.nbVar
            ratios[indexObj,indexVar] = Frac(prob.profits[indexObj,indexVar],prob.weights[indexVar])
        end
    end

    # compute the lambdas
    lambdas = Vector{LambdaChange}()
    for i in 1:prob.nbVar
        for j in (i+1):prob.nbVar
            if i == 19 && j == 38
                verbose = true
            else
                verbose = false
            end
            denominator = ratios[1,i] - ratios[2,i] - ratios[1,j] + ratios[2,j]
            if !isZero(denominator)
                numerator = ratios[2,j] - ratios[2,i]
                lambda = numerator/denominator
                if isBetweenZeroAndOne(lambda)
                    push!(lambdas,LambdaChange((i,j),lambda))
                end
            end
        end
    end

    # compute the order of the objects
    permObj = sortperm(1:prob.nbVar, lt = (x, y) -> if ratios[1, x] == ratios[1, y] isless(ratios[2, x], ratios[2, y]) else isless(ratios[1, x], ratios[1, y]) end, rev = true)

    # sort the lambdas
    lambdas = sort(lambdas, by = x->x.value, rev = true)

    # compute the "newIndicies" vector
    newIndicies = Vector{Vector{Tuple{Int,Int}}}()
    iter = 1
    nbIndiciesTotal = 0
    while iter <= length(lambdas)
        lambdaValue = lambdas[iter].value
        indicies = [lambdas[iter].index]
        differentValue = false
        if iter < length(lambdas)
            # check if the next lambda has the same value
            iter2 = iter + 1
            while iter2 <= length(lambdas) && !differentValue
                lambdaValueBuffer = lambdas[iter2].value
                if lambdaValue == lambdaValueBuffer
                    push!(indicies,lambdas[iter2].index)
                    nbIndiciesTotal += 1
                    iter += 1
                else
                    differentValue = true
                end
                iter2 += 1
            end
        end
        push!(newIndicies,indicies)
        nbIndiciesTotal += 1
        iter += 1
    end

    end # TimerOutputs

    CONFIG.debug && DEBUG_lambdas(lambdas)
    CONFIG.debug && DEBUG_newIndices(nbIndiciesTotal,length(lambdas))

    return LinearSpeedUp(newIndicies,permObj)
end

"""
Returns the relaxed UB associated with the given problem `prob` and `assignment`.

`LSU` is the `LinearSpeedUp` object provided by the initialisation.
"""
function parametricMethod(prob::BiOKP, assignment::Assignment, LSU::LinearSpeedUp)

    if assignment.lastAssigned == prob.nbVar # no unassigned vars, we have only one solution

        @timeit to "parametricMethod" begin

        x = assignment.assign[1:prob.nbVar]
        y = assignment.profit[1:2]

        sol = Sol(x,y,assignment.weight,true)

        CONFIG.debug && DEBUG_feasibleSolution(prob,sol)

        end # TimerOutput

        return UpperBound(cons(sol,nil(Sol)))

    else # still unassigned vars, we have more than one solution

        @timeit to "parametricMethod" begin

        # create the new permutation vector
        permObj = permWithoutAssignedVars(prob,LSU.permObj,assignment)

        # construct the relaxed solution associated
        sol, brokenObject, isBinary = constructRelaxedSolutionFromPermAndAssignment(prob,assignment,permObj)

        UB = UpperBound()

        UB.sols = cons(sol,nil(Sol))
        # go through each permutation couple
        for iter in 1:length(LSU.indicies)

            if length(LSU.indicies[iter]) > 1 # we need to do several permutations
                # get the permutations
                indicies = LSU.indicies[iter]
                # go through each permutation
                for iterIndex in 1:length(indicies)
                    # store the permutation couple
                    (i,j) = indicies[iterIndex]
                    # we only accept the change if i,j target unassigned variables
                    if i > assignment.lastAssigned && j > assignment.lastAssigned
                        # generate invert permutation
                        invPermObj = invertPerm(permObj,assignment.lastAssigned+1,prob.nbVar)
                        # swap i,j indicies if necessary
                        if invPermObj[j-assignment.lastAssigned] < invPermObj[i-assignment.lastAssigned]
                            i,j = j,i
                        end
                        # swap the items inside the permutation vector
                        permObj[invPermObj[i-assignment.lastAssigned]], permObj[invPermObj[j-assignment.lastAssigned]] = permObj[invPermObj[j-assignment.lastAssigned]], permObj[invPermObj[i-assignment.lastAssigned]]
                    end
                end
                # generate the new solution according to the new permutation vector
                sol,brokenObject,isBinary = constructRelaxedSolutionFromPermAndAssignment(prob,assignment,permObj)
                # if the new solution is not the same as the previous one, we add it to the list
                if sol.y != UB.sols.head.y
                    UB.sols = cons(sol,UB.sols)
                end

            else # only one permutation
                # get the permutation couple
                (i,j) = LSU.indicies[iter][1]
                # we only accept the change if i,j target unassigned variables
                if i > assignment.lastAssigned && j > assignment.lastAssigned
                    # generate invert permutation
                    invPermObj = invertPerm(permObj,assignment.lastAssigned+1,prob.nbVar)
                    # cases where the solution doesn't change (switchs are below broken object or after)
                    ijInBag = false
                    if isBinary && brokenObject == 0 # no object is broken
                        ijInBag = true
                    end
                    if !ijInBag
                        # check if the two items are before the broken item
                        ijBeforeBrokenObject = invPermObj[i-assignment.lastAssigned] < invPermObj[brokenObject-assignment.lastAssigned] && invPermObj[j-assignment.lastAssigned] < invPermObj[brokenObject-assignment.lastAssigned]
                        # check if the two iterms are after the broken item
                        if isBinary
                            ijAfterBrokenObject = invPermObj[i-assignment.lastAssigned] >= invPermObj[brokenObject-assignment.lastAssigned] && invPermObj[j-assignment.lastAssigned] >= invPermObj[brokenObject-assignment.lastAssigned]
                        else
                            ijAfterBrokenObject = invPermObj[i-assignment.lastAssigned] > invPermObj[brokenObject-assignment.lastAssigned] && invPermObj[j-assignment.lastAssigned] > invPermObj[brokenObject-assignment.lastAssigned]
                        end
                    end
                    if ijInBag || ijBeforeBrokenObject || ijAfterBrokenObject # sol is still valid, we do nothing
                        
                        # make sure that the object i is before than j in the order
                        if invPermObj[j-assignment.lastAssigned] < invPermObj[i-assignment.lastAssigned]
                            i,j = j,i
                        end
                        # switch the objects
                        permObj[invPermObj[i-assignment.lastAssigned]], permObj[invPermObj[j-assignment.lastAssigned]] = permObj[invPermObj[j-assignment.lastAssigned]], permObj[invPermObj[i-assignment.lastAssigned]]
                        # update the broken item
                        if brokenObject == i
                            brokenObject = j
                        elseif brokenObject == j
                            brokenObject = i
                        end

                    else # sol isn't valid anymore, we construct and add new sol to the list
                        
                        # generate invert permutation vector
                        invPermObj = invertPerm(permObj,assignment.lastAssigned+1,prob.nbVar)
                        # make sure that the object i is before than j in the order
                        if invPermObj[j-assignment.lastAssigned] < invPermObj[i-assignment.lastAssigned]
                            i,j = j,i
                        end
                        # switch the objects
                        permObj[invPermObj[i-assignment.lastAssigned]], permObj[invPermObj[j-assignment.lastAssigned]] = permObj[invPermObj[j-assignment.lastAssigned]], permObj[invPermObj[i-assignment.lastAssigned]]
                        # construct the solution
                        sol,brokenObject,isBinary = constructRelaxedSolutionFromPermAndAssignment(prob,assignment,permObj)

                        UB.sols = cons(sol, UB.sols)

                    end
                end
            end
        end

        end # TimerOutput

        CONFIG.debug && DEBUG_UB(prob, UB)

        return UB
    end
end

#-------------------------------------------------------------------#
#----------------------- PRIVATE FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

"""
Returns the invert permutation vector `rev` of the permutation vector `perm`.

Used to access in constant time to an element in the `perm` list.

# Example

```jldoctest
julia> invertPerm([3,1,2],1,3)
6-element Array{Int64,1}:
 2
 3
 1
```
"""
function invertPerm(perm::Vector{Int}, minValue::Int, maxValue::Int)
	rev = Vector{Int}(undef, maxValue - minValue + 1)
	for iter in 1:(maxValue - minValue + 1)
		rev[perm[iter]-minValue+1] = iter
	end
	return rev
end

"""
Returns a relaxed solution, filling the bag according to the given permutation vector `permObj`and `assignment`.

`assignment` is the assignment (stating which variables are fixed).

`permObj` is the permutation vector (order didacting which variables go in the bag first).
"""
function constructRelaxedSolutionFromPermAndAssignment(prob::BiOKP, assignment::Assignment, permObj::Vector{Int})

    CONFIG.debug && @assert length(permObj) == prob.nbVar - assignment.lastAssigned

    @timeit to "constructRelaxed" begin
    
    # create the partial solution associated with the assignment
    weight = assignment.weight
    y = copy(assignment.profit)
    x = zeros(Float64,prob.nbVar)
    for i in 1:assignment.lastAssigned
        x[i] = assignment.assign[i]
    end

    i = 1
    # fill the bag
    while  i <= prob.nbVar-assignment.lastAssigned && prob.maxWeight - weight >= prob.weights[permObj[i]] # while we can put the next object
        x[permObj[i]] = 1 # we add the object
        weight += prob.weights[permObj[i]] # we substract the weight
        y += [prob.profits[1,permObj[i]],prob.profits[2,permObj[i]]]
        i += 1
    end

    # put the broken object
    brokenObjectIndex = i
    isBinary = (weight == prob.maxWeight) || brokenObjectIndex == (length(permObj)+1)
    if !isBinary
        @assert brokenObjectIndex <= (length(permObj) + 1) "Out of range"
        x[permObj[brokenObjectIndex]] = (prob.maxWeight - weight)/prob.weights[permObj[brokenObjectIndex]]
        y[1] += (prob.maxWeight - weight)/prob.weights[permObj[brokenObjectIndex]]*prob.profits[1,permObj[brokenObjectIndex]]
        y[2] += (prob.maxWeight - weight)/prob.weights[permObj[brokenObjectIndex]]*prob.profits[2,permObj[brokenObjectIndex]]
        weight = prob.maxWeight
    end

    sol = Sol(x,y,weight,isBinary)

    end # TimerOutput

    CONFIG.debug && DEBUG_feasibleSolution(prob, sol)

    return sol, brokenObjectIndex <= (length(permObj)) ? permObj[brokenObjectIndex] : 0, isBinary
end

"""
Return a permutation vector `returnedPermObj` taking into account only the unassigned variables.

`permObj` is a full permutation vector.

`assignment` is the current assignment.
"""
function permWithoutAssignedVars(prob::BiOKP, permObj::Vector{Int}, assignment::Assignment)

    @timeit to "permWithoutAssignedVars" begin

    # permutation vector but without the assigned variables
    returnedPermObj = Vector{Int}(undef,prob.nbVar-assignment.lastAssigned)

    iter = 1
    nbVarsPlaced = 0
    while iter <= prob.nbVar && nbVarsPlaced < prob.nbVar - assignment.lastAssigned
        # if the var is unassigned, we put it inside the new permutation vector
        if permObj[iter] > assignment.lastAssigned
            nbVarsPlaced += 1
            returnedPermObj[nbVarsPlaced] = permObj[iter]
        end
        iter += 1
    end

    end # TimerOutputs

    return returnedPermObj
end