
""" Upper bound computation (prefix is UB for Upper Bound) """

# All functions are public :

"""
Returns an object `parentToChild` : the already known points that won't need to be computed again.

Points we already know are the ones for which the variable `index` is already set to `value`.

Returns parentToChild
"""
function UB_computeParentToChild(prob::BiOKP, UB::UpperBound, index::Int, value::Int)

    CONFIG.debug && DEBUG_UB(prob, UB)

	@timeit to "UB_computeParentToChild" begin

    # known solutions that go from the parent to the child node
	knownSols = nil(Sol)
	knownSolsLastElt = nil(Sol)

    # indicies of the known solutions inside the parent UB.sols
    previousIndicies = nil(Int)
    previousIndiciesLastElt = nil(Int)

    sols = UB.sols

    lenLB = length(sols)

    # attributes of the ParentToChild object : are lexicographic sols known by the child ?
    leftKnown = false
    rightKnown = false

    # right sol to memorize
    rightSol = Sol()

    iter = 1
    atLeastOneSolKnown = false # length(knownSols) >= 1
	while sols != nil(Sol)

		if !atLeastOneSolKnown && sols.head.x[index] == value

            # length(knownSols) >= 1 so atLeastOneSolKnown becomes true
			atLeastOneSolKnown = true

            # append the known solution to knownSols
			knownSols = cons(sols.head, nil(Sol))
			knownSolsLastElt = knownSols

            # append the index of the known solution to previousIndicies
            previousIndicies = cons(iter, nil(Int))
			previousIndiciesLastElt = previousIndicies

            rightSol = sols.head

            # update leftKnown and rightKnown booleans
            if iter == 1 # the left sol is a known solution for the child
                leftKnown = true
            elseif iter == lenLB
                rightKnown = true
            end

		elseif sols.head.x[index] == value

            # append the known solution to knownSols
			knownSolsLastElt.tail = cons(sols.head, nil(Sol))
			knownSolsLastElt = knownSolsLastElt.tail

            # append the index of the known solution to previousIndicies
            previousIndiciesLastElt.tail = cons(iter, nil(Int))
			previousIndiciesLastElt = previousIndiciesLastElt.tail

            rightSol = sols.head

            # update leftKnown and rightKnown booleans
            if iter == 1 # the left sol is a known solution for the child
                leftKnown = true
            elseif iter == lenLB
                rightKnown = true
            end

		end

		sols = sols.tail
        iter += 1
	end

    parentToChild = ParentToChild(knownSols,previousIndicies,rightSol,leftKnown,rightKnown)

    end # TimerOutput

    CONFIG.debug && DEBUG_parentToChild(prob, parentToChild, index=index, value=value)

	return parentToChild
end

"""
Computes the segments associated with the given `UB.sols` as a DualSet
"""
function UB_computeSegments!(prob::BiOKP, UB::UpperBound)

    CONFIG.debug && DEBUG_UB(prob, UB)

    @timeit to "UB_computeSegments!" begin

    sols = UB.sols

	nbConstraints = length(sols) + 1

	A = zeros(Float64, nbConstraints, 2)
	b = zeros(Float64, nbConstraints)

	max1 = sols.head.y[1]
	max2 = sols.head.y[2]

	iter = 1

	while sols.tail != nil(Sol)

        solL = sols.head
        solR = sols.tail.head

		A[iter, 1] = solL.y[2] - solR.y[2]
		A[iter, 2] = solR.y[1] - solL.y[1]

		b[iter] = (solL.y[2] - solR.y[2]) * solR.y[1] + (solR.y[1] - solL.y[1]) * solR.y[2]

		if solL.y[1] > max1
			max1 = solL.y[1]
		end
		if solR.y[1] > max1
			max1 = solR.y[1]
		end
		if solL.y[2] > max2
			max2 = solL.y[2]
		end
		if solR.y[2] > max2
			max2 = solR.y[2]
		end

		sols = sols.tail
		iter += 1
	end

    # add the constraint extreme horizontal and vertical to finish the polytope

    # vertical (on the right)
	A[nbConstraints-1, 1] = 1
	A[nbConstraints-1, 2] = 0

	b[nbConstraints-1] = max1

    # horizontal (on the left)
	A[nbConstraints, 1] = 0
	A[nbConstraints, 2] = 1

	b[nbConstraints] = max2

    UB_setSegments!(UB, DualSet(A, b)) # update segments

    end # TimerOutput

    CONFIG.debug && DEBUG_UB(prob, UB)

end