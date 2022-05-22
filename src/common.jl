
""" Common operations made in the code """

# All functions are public

"""
Add the given solution `sol` to the current `front`.

UNSAFE : `sol` is assumed to be non dominated in `front` and to not dominate any solution. Please use this function only when `sol` verifies these conditions (in a dichotomic search for example). Otherwise the code will throw an error.

The solution `sol` is inserted in a way that the `front` remains sorted and keeps unicity.
"""
function strictInsertInFront!(prob::BiOKP, front::LinkedList{Sol}, sol::Sol)

    CONFIG.debug && DEBUG_front(prob, front)
    CONFIG.debug && DEBUG_feasibleSolution(prob, sol)

    @timeit to "strictInsertInFront!" begin

	solHasBeedInserted = false
    # going through the linked list to insert the solution at the right place
	while front.tail != nil(Sol) && !solHasBeedInserted
		if front.tail.head.y[1] > sol.y[1]
			front.tail = cons(sol, front.tail)
			solHasBeedInserted = true
		end
		front = front.tail
	end
    # handling a possible bad case scenario
	@assert solHasBeedInserted "The solution couldn't be inserted. This not usual.\n$front\n$(sol.y)"

    end # TimerOutputs

    CONFIG.debug && DEBUG_front(prob, front)
end

"""
Merge the given two lists of nadir points : `nadirs1` and `nadirs2`.
"""
function mergeNadirs(prob::BiOKP, nadirs1::Vector{PairOfSolution}, nadirs2::Vector{PairOfSolution})

    CONFIG.debug && DEBUG_nadirs(prob, nadirs1)
    CONFIG.debug && DEBUG_nadirs(prob, nadirs2)

    @timeit to "mergeNadirs" begin

	lengthA = length(nadirs1)
	lengthB = length(nadirs2)

    # extreme cases : one or both lists are empty
    if lengthA == 0 && lengthB == 0 # both
        return nadirs1
    elseif lengthA == 0
        return nadirs2
    elseif lengthB == 0
        return nadirs1
    end

	iterA = 1 # iterator on list A
	iterB = 1 # iterator on list B
	iterF = 1 # iterator on list F

    # allocate the result
	finalList = Vector{PairOfSolution}(undef, lengthA + lengthB)

    # number of elements in the result
    nbTrueElts = 0

    # fill the elements in order until one iterator reaches the end of its list
	while iterA <= lengthA && iterB <= lengthB
		if nadirs1[iterA].solL.y[1] <= nadirs2[iterB].solL.y[1] # element A is chosen
			finalList[iterF] = nadirs1[iterA]
            nbTrueElts += 1
			iterA += 1
            iterF += 1
		else
            # check if the nadir is a different one than the previous one
            # if so insert element B
            if iterF > 1 && finalList[iterF-1].solL.y[1] != nadirs2[iterB].solL.y[1]
                finalList[iterF] = nadirs2[iterB]
                nbTrueElts += 1
                iterF += 1
            elseif iterF == 1
                finalList[iterF] = nadirs2[iterB]
                nbTrueElts += 1
                iterF += 1
            end
            iterB += 1
		end
	end

    # finish by inserting the remaining elements
	if iterA > lengthA
        for i in iterB:lengthB
            if nbTrueElts+1 > 1 && finalList[nbTrueElts].solL.y[1] != nadirs2[iterB].solL.y[1]
                finalList[nbTrueElts+1] = nadirs2[iterB]
                nbTrueElts += 1
            end
            iterB += 1
        end
	else
		for i in iterA:lengthA
            finalList[nbTrueElts+1] = nadirs1[iterA]
            iterA += 1
            nbTrueElts += 1
        end
	end

    # keep only an array of the right length
    finalList = finalList[1:nbTrueElts]

    end # TimerOutputs

    CONFIG.debug && DEBUG_nadirs(prob, finalList)

	return finalList
end

#--------------------- SOLVERS GENERATORS ---------------------#

" Returns the solver selected for the first LB_dichoSearch "
function getSolver1OKP(firstDicho::FirstDicho)
    if firstDicho == JUMP
        return SOLVE1OKP_GLPK
    elseif firstDicho == COMBO
        return SOLVE1OKP_combo
    else
        @assert false "The method is not supported"
    end
end

" Returns the solver selected for the computation of the UB "
function getSolver1OKP(methodUB::MethodUB)
    if methodUB == EXACT_JUMP
        return SOLVE1OKP_GLPK
    elseif methodUB == EXACT_COMBO
        return SOLVE1OKP_combo
    elseif methodUB == RELAX_LIN_CLASSIC
        return SOLVE1OKP_linear
    else
        @assert false "The method is not supported"
    end
end

#------------------ USED in DICHOSEARCH and SubProb Study --------------

"""
Returns the optimal lexigraphic solutions for the given problem `prob` and `assignment`.

Giving a hydrated `parentToChild` object helps fasten the computation. In some cases we can avoid solving some scalarized problems.
"""
function getLexicoSolutions(prob::BiOKP, assignment::Assignment, frontOfParent::LinkedList{Sol}, solver1OKP::Function; parentToChild::ParentToChild = ParentToChild())

    CONFIG.debug && DEBUG_front(prob, frontOfParent)
    CONFIG.debug && DEBUG_parentToChild(prob, parentToChild)

    @timeit to "getLexicoSolutions" begin

    # Compute the first lexicographic optimal solution
    if !parentToChild.leftKnown # the left sol isn't known by parentToChild
        λ = [1, 1000.]
        leftSol = solver1OKP(prob, λ, assignment)
    else
        leftSol = frontOfParent.head
    end

    # Compute the second lexicographic optimal solution
    if !parentToChild.rightKnown # the right sol isn't known by parentToChild
        λ = [1000., 1]
        rightSol = solver1OKP(prob, λ, assignment)
    else
        rightSol = parentToChild.rightSol
    end

    # update parentToChild object

    if !parentToChild.leftKnown && parentToChild.knownSols != nil(Sol) && parentToChild.knownSols.head.y == leftSol.y
        parentToChild.leftKnown = true
    end

    if !parentToChild.rightKnown && parentToChild.knownSols != nil(Sol) && parentToChild.rightSol.y == rightSol.y
        parentToChild.rightKnown = true
    end

    end # TimerOutputs

    CONFIG.debug && DEBUG_feasibleSolution(prob, leftSol)
    CONFIG.debug && DEBUG_feasibleSolution(prob, rightSol)

    return leftSol,rightSol
end

"""
Returns the pairs of solutions `to study` during the dichotomic search and the given `front` of solutions completed if possible with lexicographic solutions.

The solutions `leftSol` and `rightSol` are the lexicographic optimal solutions for the current sub problem.

`NOTE` : this function comes from an optimization we made to fasten the dichotomic method. At first, the dichotomic search was studying all the pairs available. Now we avoid the pairs of solutions that were consecutive inside the parent node.
"""
function computesPairsToStudy(prob::BiOKP, front::LinkedList{Sol}, parentToChild::ParentToChild, leftSol::Sol, rightSol::Sol)

    CONFIG.debug && DEBUG_front(prob, front)
    CONFIG.debug && DEBUG_parentToChild(prob, parentToChild)
    CONFIG.debug && DEBUG_feasibleSolution(prob, leftSol)
    CONFIG.debug && DEBUG_feasibleSolution(prob, rightSol)

    @timeit to "computesPairsToStudy" begin

    toStudy = Vector{PairOfSolution}()

    # add the left sol if it hasn't been done yet
    if front == nil(Sol) || (!parentToChild.leftKnown && leftSol.y != front.head.y)
        front = cons(leftSol, front)
    end

    # different scenarios for different sizes of parentToChild.knownSols
    if parentToChild.knownSols == nil(Sol) # no known solutions
        # add rightSol at the end of the front
        front.tail = cons(rightSol,nil(Sol))
        # add the pair left,right
        push!(toStudy, PairOfSolution(front.head, front.tail.head))
    elseif parentToChild.knownSols.tail == nil(Sol) # one known solution
        if parentToChild.leftKnown || parentToChild.rightKnown # this known solution is the left one or right one
            # add rightSol at the end of the front
            front.tail = cons(rightSol,nil(Sol))
            # add the pair left,right
            push!(toStudy, PairOfSolution(front.head, front.tail.head))
        else # this known solution is in the middle
            # add rightSol at the end of the front
            front.tail.tail = cons(rightSol,nil(Sol))
            # add the pair left,middle
            push!(toStudy, PairOfSolution(front.head, front.tail.head))
            # add the pair middle,right
            push!(toStudy, PairOfSolution(front.tail.head, front.tail.tail.head))
        end
    else # more than one solution known

        # rightSol wasn't in the knownSols list, thus isn't in the current front
        if !parentToChild.rightKnown
            ll = front
            while ll.tail != nil(Sol)
                ll = ll.tail
            end
            # add rightSol at the end of the front
            ll.tail = cons(rightSol,nil(Sol))
        end

        # leftSol wasn't in the knownSols list, in any case we'll have to study the pair (left,firstKnownSol)
        if !parentToChild.leftKnown
            push!(toStudy, PairOfSolution(front.head, front.tail.head))
        end

        # adding the remaining pairs when necessary
        prevInd = parentToChild.previousIndicies
        ll = front
        while prevInd.tail != nil(Sol)
            # if the known solutions weren't consecutive for the parent, we have to study the pair
            if prevInd.head + 1 != prevInd.tail.head
                # leftKnown value shifts the front, we take this into account here to study the right pair of solutions
                if parentToChild.leftKnown # the list is the same length as prevInd
                    push!(toStudy, PairOfSolution(ll.head,ll.tail.head))
                else # the list is shifted from one to the right
                    push!(toStudy, PairOfSolution(ll.tail.head,ll.tail.tail.head))
                end
            end

            # if the right solution wasn't in the parent, we always need to add the pair ("sol before",rightSol). But we need to do this at the end of the while loop, that's why we make sure that "prevInd.tail.tail == nil(Sol)".
            if !parentToChild.rightKnown && prevInd.tail.tail == nil(Sol)
                # again we take care of the shift caused by leftKnown boolean
                if parentToChild.leftKnown # no shift
                    push!(toStudy, PairOfSolution(ll.tail.head,ll.tail.tail.head))
                else # shift
                    push!(toStudy, PairOfSolution(ll.tail.tail.head,ll.tail.tail.tail.head))
                end
            end

            prevInd = prevInd.tail
            ll = ll.tail
        end
    end

    end # TimerOutput

    CONFIG.debug && DEBUG_front(prob, front)
    CONFIG.debug && DEBUG_pairsOfSolutions(prob, toStudy)

    return toStudy, front
end