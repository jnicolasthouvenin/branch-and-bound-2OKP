
""" Functions indicating if a subproblem can be fathomed and giving information on the dominated nadir points of the lower bound and upper Bound solutions """

# Public functions :
# - SPS_partially_implicit : study sub problem using partial dichotomy

#-------------------------------------------------------------------#
#------------------------ PUBLIC FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

"""
Optimized pruning test (for dichoSearch type UB computation)

Computes the UB of the subProblem and stops as soon as the sub problem cannot be fathomed.
"""
function SPS_partially_implicit(prob::BiOKP, assignment::Assignment, nadirsToStudy::Vector{PairOfSolution}; parentToChild = ParentToChild())

    CONFIG.debug && DEBUG_nadirs(prob,nadirsToStudy)
    CONFIG.debug && DEBUG_parentToChild(prob,parentToChild)

    @timeit to "subProblemStudy" begin

    UB = UpperBound() # init

    # get initial solutions
    UB.sols = copy(parentToChild.knownSols)

    # compute lexicographic optimal solutions
    leftSol, rightSol = getLexicoSolutions(prob, assignment, UB.sols, getSolver1OKP(COMPONENTS.methodUB), parentToChild = parentToChild)

    # In the event that the two solutions found are identical, this solution is optimal and thus added to the lower bound
    if !differentFloatSols(leftSol,rightSol)
        if UB.sols == nil(Sol) # if the UB.sols contains already one solution, we don't need to add this one because it is the same
            UB.sols = cons(leftSol, UB.sols)
        end
        # we don't have any pairs to study

        if leftSol.isBinary
            CONFIG.debug && DEBUG_UB(prob, UB)
            return OPTIMALITY, nadirsToStudy, Vector{PairOfSolution}(), UB
        else
            CONFIG.debug && DEBUG_UB(prob, UB)
            return NONE, nadirsToStudy, Vector{PairOfSolution}(), UB
        end

    else # the two lexicographic solutions are different
        
        # get the list of pairs of solutions to study
        toStudy, UB.sols = computesPairsToStudy(prob, UB.sols, parentToChild, leftSol, rightSol)

        # add the constraint extreme horizontal and vertical of the polytope
        A = Matrix{Float64}(undef,2,2)
        b = Vector{Float64}(undef,2)

        # vertical (on the right)
        A[1, 1] = 1
        A[1, 2] = 0

        b[1] = rightSol.y[1]

        # horizontal (on the left)
        A[2, 1] = 0
        A[2, 2] = 1

        b[2] = leftSol.y[2]

        subDominated, nadirsToStudy, nonDomNadirs = PRUNING_compareCstrToNadirs(prob, DualSet(A,b), nadirsToStudy)

        # while some pairs of solutions are found by the dichotomic search, we keep going
        while toStudy != [] && !subDominated

            currPair = pop!(toStudy)

            leftSol = currPair.solL
            rightSol = currPair.solR

            λ = [leftSol.y[2] - rightSol.y[2], rightSol.y[1] - leftSol.y[1]]

            solver1OKP = getSolver1OKP(COMPONENTS.methodUB)
            midSol     = solver1OKP(prob, λ, assignment)

            # if the solution is beyond the pair studied, it's added to the UB
            if round(sum(λ .* midSol.y), digits = 7) - round(sum(λ .* rightSol.y), digits = 7) > 0.00001
                # add the new solution found
                strictInsertInFront!(prob, UB.sols, midSol)
                # add the two newly found pairs
                push!(toStudy, PairOfSolution(leftSol, midSol))
                push!(toStudy, PairOfSolution(midSol, rightSol))
            end

            # we construct the constraint associated with the new solution
            b = (leftSol.y[2] - rightSol.y[2]) * midSol.y[1] + (rightSol.y[1] - leftSol.y[1]) * midSol.y[2]

            A = Matrix{Float64}(undef,1,2)
            A[1,1] = λ[1]
            A[1,2] = λ[2]

            subDominated, nadirsToStudy, nonDomNadirs = PRUNING_compareCstrToNadirs(prob, DualSet(A,[b]),nadirsToStudy, nonDomNadirs = nonDomNadirs)
        end

        if subDominated

            CONFIG.debug && DEBUG_nadirs(prob,nadirsToStudy)
            CONFIG.debug && DEBUG_UB(prob, UB)

            return DOMINANCE, Vector{PairOfSolution}(), nadirsToStudy, UB

        else

            CONFIG.debug && DEBUG_nadirs(prob, nadirsToStudy)
            CONFIG.debug && DEBUG_nadirs(prob, nonDomNadirs)
            CONFIG.debug && DEBUG_UB(prob, UB)

            return NONE, nadirsToStudy, nonDomNadirs, UB

        end
    end

    end # TimerOutput
end

#-------------------------------------------------------------------#
#----------------------- PRIVATE FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

"""
Indicates if the subproblem is dominated, as well as a partition of nadir points (dominated vs non dominated)

[ATTENTION] For nadir points, "dominated" doesn't mean "discarded". It is the OPPOSITE ! If the subproblem isn't dominated, we will keep studying the dominated nadirs, not the nondominated ones. In the rest of the code, dominated nadirs are therefore called "nadirsToStudy".

Returns isSubDominated, dominatedNadirs, nonDomNadirs
"""
function PRUNING_compareCstrToNadirs(prob::BiOKP, constraint::DualSet, nadirsToStudy::Vector{PairOfSolution}; nonDomNadirs::Vector{PairOfSolution} = Vector{PairOfSolution}(), regionLimit::Vector{Float64} = [-1.,-1.])

    CONFIG.debug && DEBUG_nadirs(prob,nadirsToStudy)
    CONFIG.debug && DEBUG_nadirs(prob,nonDomNadirs)

    @timeit to "PRUNING_compareCstrToNadirs" begin
    
    dominatedNadirs = Vector{PairOfSolution}()

    regionLimit == [-1,-1] ? checkLimits = false : checkLimits = true

    A = constraint.A
    b = constraint.b

    isSubDominated = true
    iterNadirs = 0
    # for each nadir point to study
    while iterNadirs < length(nadirsToStudy)

        iterNadirs += 1

        pairNadir = nadirsToStudy[iterNadirs]

        if prob.isInteger && COMPONENTS.nadirsShift
		    nadir = [pairNadir.solL.y[1] + 1, pairNadir.solR.y[2] + 1] # - 0.0...1 to avoid bad roundings
        else
            nadir = [pairNadir.solL.y[1], pairNadir.solR.y[2]]
        end

        # if the nadir is out of the bounds
        if checkLimits && (nadir[1] > regionLimit[1] || nadir[2] > regionLimit[2])
            # the nadir is non dominated
            push!(nonDomNadirs,pairNadir)

        else

            nadirA = A * nadir

            isNadirUnderUB = true
            iterConstraint = 1
            while iterConstraint <= length(b) && isNadirUnderUB
                # compute isNadirUnderUB : verify that the nadir point is lower that the upper bound on all objectives
                if prob.isInteger && COMPONENTS.nadirsShift
                    isNadirUnderUB = isNadirUnderUB && nadirA[iterConstraint] <= b[iterConstraint]
                else
                    isNadirUnderUB = isNadirUnderUB && nadirA[iterConstraint] < b[iterConstraint]
                end
                iterConstraint += 1
            end

            if isNadirUnderUB # the nadir is dominated but the UB
                push!(dominatedNadirs,pairNadir)
                isSubDominated = false
            else
                push!(nonDomNadirs,pairNadir)
            end
        end
	end

    while iterNadirs < length(nadirsToStudy)
        iterNadirs += 1
        push!(dominatedNadirs,nadirsToStudy[iterNadirs])
    end

    sort!(nonDomNadirs, lt = (x, y) -> isless(x.solL.y[1],y.solL.y[1]))

    end # TimerOutput

    CONFIG.debug && DEBUG_nadirs(prob, dominatedNadirs)
    CONFIG.debug && DEBUG_nadirs(prob, nonDomNadirs)

    return isSubDominated, dominatedNadirs, nonDomNadirs
end