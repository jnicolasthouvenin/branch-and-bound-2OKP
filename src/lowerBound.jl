
""" Working with the LowerBound (prefix is LB for Lower Bound) """

# Public functions :
# - LB_dichoSearch : computes X_SE1
# - LB_update! : add potential new points to the given LB
# - LB_getNadirs : computes the local nadir points of a given LB

#-------------------------------------------------------------------#
#------------------------ PUBLIC FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

"""
Returns the set LowerBound object (`X_SE1`) associated with the given problem `prob` and `assignment`, computed by the dichotomic method.

Two solvers are available : `Combo` (EXACT_COMBO) and `GLPK` (EXACT_JUMP).

To speed up the search, a parameter `parentToChild` containing a list of already known solutions can be provided.
"""
function LB_dichoSearch(prob::BiOKP, assignment::Assignment; parentToChild = ParentToChild())

    CONFIG.debug && DEBUG_parentToChild(prob, parentToChild)

    LB = LowerBound() # init empty lower bound

    LB.sols = copy(parentToChild.knownSols)

    # compute left and right solutions
    leftSol, rightSol = getLexicoSolutions(prob, assignment, LB.sols, getSolver1OKP(COMPONENTS.firstDicho), parentToChild = parentToChild)

    # In the event that the two solutions found are identical, this solution is optimal and thus added to the lower bound
    if !differentFloatSols(leftSol,rightSol)
        if LB.sols == nil(Sol) # if the LB contains already one solution, we don't need to add this one because it is the same
            LB.sols = cons(leftSol, LB.sols)
        end
    else # the two lexicographic solutions are different

        # get the list of pairs of solutions to study
        toStudy, LB.sols = computesPairsToStudy(prob, LB.sols, parentToChild, leftSol, rightSol)

		# while some pairs of solutions are found by the dichotomic search, we keep going
        while toStudy != []

			currPair = pop!(toStudy)

			leftSol = currPair.solL
			rightSol = currPair.solR

	        λ = [leftSol.y[2] - rightSol.y[2], rightSol.y[1] - leftSol.y[1]]

            solver1OKP = getSolver1OKP(COMPONENTS.firstDicho)
            midSol = solver1OKP(prob, λ, assignment)

			# if the solution is beyond the pair studied, it's added to the LB
			if round(sum(λ .* midSol.y), digits = 7) - round(sum(λ .* rightSol.y), digits = 7) > 0.00001
                # add the new solution found
				strictInsertInFront!(prob, LB.sols, midSol)
                # add the two newly found pairs
				push!(toStudy, PairOfSolution(leftSol, midSol))
				push!(toStudy, PairOfSolution(midSol, rightSol))
            end
        end
    end

    CONFIG.debug && DEBUG_LB(prob, LB)

	return LB

end

"""
Update the given `LB` and returns the new vector of nadir points to study `newNadirsToStudy` according to the given `UBSols` and `dominatedNadirs`.

`dominatedNadirs` are the nadir points to study.
"""
function LB_update!(prob::BiOKP, LB::LowerBound, dominatedNadirs::Vector{PairOfSolution}, UB::UpperBound)

    CONFIG.debug && DEBUG_LB(prob, LB)
    CONFIG.debug && DEBUG_UB(prob, UB)
    CONFIG.debug && DEBUG_nadirs(prob, dominatedNadirs)

    for sol in UB.sols
        # check if the solution found in the upper bound is binary (thus feasible)
        if !(COMPONENTS.methodUB == RELAX_LIN_CLASSIC) || sol.isBinary

            LB_removeAllDominatedSols!(prob, LB, sol)

            # insert sol at the right place (if the solution is not already dominated by the current LB)
            inserted = LB_insert!(prob, LB, sol)

        end
    end

    newNadirsToStudy = LB_getNadirsWithNadirsToStudy(prob, LB, dominatedNadirs)

    CONFIG.debug && DEBUG_LB(prob, LB)
    CONFIG.debug && DEBUG_nadirs(prob, newNadirsToStudy)

    return newNadirsToStudy
end

"""
Returns the nadir points `nadirs` associated with the given `LB`.
"""
function LB_getNadirs(prob::BiOKP, LB::LowerBound)

    CONFIG.debug && DEBUG_LB(prob, LB)

    sols = LB.sols
    nadirs = Vector{PairOfSolution}(undef, length(sols) - 1)

    iterNadirs = 1
    # go through sols and simply create each nadir point
    while sols.tail != nil(Sol)
    	nadirs[iterNadirs] = PairOfSolution(sols.head, sols.tail.head)
    	sols = sols.tail
    	iterNadirs += 1
    end

    CONFIG.debug && DEBUG_nadirs(prob, nadirs)

    return nadirs
end

#-------------------------------------------------------------------#
#----------------------- PRIVATE FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

"""
Inserts in place (if possible) the given solution `sol` into the given `LB`.

The solution `sol` is assumed to `NOT DOMINATE` solutions inside `LB.sols`.
"""
function LB_insert!(prob::BiOKP, LB::LowerBound, sol::Sol)

    CONFIG.debug && DEBUG_LB(prob, LB)
    CONFIG.debug && DEBUG_feasibleBinarySolution(prob, sol)
    CONFIG.debug && DEBUG_isntDominatingAny(sol, LB.sols)

    ll = LB.sols
    while ll != nil(Sol)
        # sol is dominated by the head of LB.sols or is the same, we shouldn't insert it
        if dominate(ll.head,sol) || (prob.isInteger && sol.y == ll.head.y) || (!prob.isInteger && !differentFloatSols(sol,ll.head))
            CONFIG.debug && DEBUG_LB(prob, LB)
            return false

        # we have to insert before the head
        elseif sol.y[1] < ll.head.y[1]
            oldSol = ll.head
            ll.head = sol
            CONFIG.debug && DEBUG_LB(prob, LB)
            return LB_insert!(prob, ll, oldSol)

        else

            if ll.tail != nil(Sol) # we can keep going

                # sol is dominated by the following solution or is the same, we shouldn't insert it
                if dominate(ll.tail.head,sol) || (prob.isInteger && sol.y == ll.tail.head.y) || (!prob.isInteger && !differentFloatSols(sol,ll.tail.head))
                    CONFIG.debug && DEBUG_LB(prob, LB)
                    return false

                # we have to insert sol between ll and ll.tail
                elseif sol.y[1] < ll.tail.head.y[1] # we have to insert sol between ll and ll.tail
                    oldTail = ll.tail
                    ll.tail = cons(sol,oldTail)
                    CONFIG.debug && DEBUG_LB(prob, LB)
                    return true
                end

            else # we are at the end of the list, we have to insert sol after ll
                ll.tail = cons(sol,nil(Sol))
                CONFIG.debug && DEBUG_LB(prob, LB)
                return true
            end

        end

        ll = ll.tail
    end

    CONFIG.debug && DEBUG_LB(prob, LB)

    return false
end

"""
Remove in place from the given `LB`, solutions that are dominated by the given solution `sol`, and returns wether or not a solution from `LB` has been removed.
"""
function LB_removeAllDominatedSols!(prob::BiOKP, LB::LowerBound, sol::Sol)

    CONFIG.debug && DEBUG_LB(prob, LB)
    CONFIG.debug && DEBUG_feasibleBinarySolution(prob, sol)

    function LB_removeOneDominatedSol!(prob::BiOKP, front::LinkedList{Sol}, sol::Sol)

        ll = front

        removed = false
        if dominate(sol,ll.head)
            if ll.tail != nil(Sol)
                if ll.tail.tail != nil(Sol)
                    ll.head = ll.tail.head
                    ll.tail = ll.tail.tail
                else
                    ll.head = ll.tail.head
                    ll.tail = nil(Sol)
                end
                removed = true
            end
        end

        if ll != nil(Sol) && ll.tail != nil(Sol) && ll.tail.tail == nil(Sol)
            if dominate(sol,ll.tail.head)
                removed = true
                ll.tail = nil(Sol)
            end
        end

        CONFIG.debug && DEBUG_front(prob, front)

        return removed

    end

    # delete in LB the solutions dominated by sol
    front = LB.sols
    while front != nil(Sol)
        removed = LB_removeOneDominatedSol!(prob, front, sol)
        if !removed
            front = front.tail
        end
    end

    CONFIG.debug && DEBUG_LB(prob, LB)

end

"""
Returns the nadir points `nadirs` associated with the given `LB`, but this time we care only about nadirs that are within the `nadirsToStudy` list.
"""
function LB_getNadirsWithNadirsToStudy(prob::BiOKP, LB::LowerBound, nadirsToStudy::Vector{PairOfSolution})

    CONFIG.debug && DEBUG_LB(prob, LB)
    CONFIG.debug && DEBUG_nadirs(prob, nadirsToStudy)

    # if the nadirsToStudy list is empty, just call the default function LB_getNadirs
    if length(nadirsToStudy) == 0
        return LB_getNadirs(prob, LB)
    end

    nadirs = Vector{PairOfSolution}()

    sols = LB.sols

    iterNadirs = 1
    # go through the sols and push certain nadirs
    while sols.tail != nil(Sol) && iterNadirs <= length(nadirsToStudy)

        # store the two current solutions from the sols
        solLBLeft = sols.head
        solLBRight = sols.tail.head
        # store the coordinates of the nadir associated with the two solutions
        z1,z2 = solLBLeft.y[1],solLBRight.y[2]

        # check if the nadir is to be stored
        if z1 >= nadirsToStudy[iterNadirs].solL.y[1] && z2 >= nadirsToStudy[iterNadirs].solR.y[2]
            push!(nadirs,PairOfSolution(solLBLeft,solLBRight))
            sols = sols.tail

        # check if the nadir won't even be stored
        elseif z1 < nadirsToStudy[iterNadirs].solL.y[1] && z2 > nadirsToStudy[iterNadirs].solR.y[2]
            sols = sols.tail

        # check if the nadir will eventually be stored
        elseif z1 > nadirsToStudy[iterNadirs].solL.y[1] && z2 < nadirsToStudy[iterNadirs].solR.y[2]
            iterNadirs += 1
        end
    end

    CONFIG.debug && DEBUG_nadirs(prob, nadirs)

    return nadirs
end