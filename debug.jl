
""" Assertions on the variables used in the program. Helps locating errors right away """

# All functions are public

# All debug functions run if the global parameter CONFIG.debug is set to true.
# These functions are expensive. When launching benchmarks, make sure to set
# CONFIG.debug to false.

# The file is divided in two sections :
#   - core : the debug functions used in the core source files
#   - parametric relaxation : the debug functions specific to the parametric relax (optional)

#------------------------- CORE -------------------------#

" `Assertion` : `correct`, `feasible`, `integer`, `sorted z1 increasing`, `unicity`, `non dominated solutions`"
function DEBUG_LB(prob::BiOKP, LB::LowerBound)

    # Change type of LB.sols
    if typeof(LB.sols) != Vector{Sol}
        # convert into array (easier to work with indices)
        sols = Vector{Sol}()
        for sol in LB.sols
            push!(sols, sol)
        end
    else
        sols = LB.sols
    end
    
    for i in 1:length(sols)

        # feasible
        DEBUG_feasibleBinarySolution(prob, sols[i])

        for j in (i+1):length(sols)

            # unicity
            if prob.isInteger
                @assert sols[i].y != sols[j].y "The LB contains duplicates\n$sols"
            else
                @assert differentFloatSols(sols[i],sols[j]) "The LB contains duplicates ($(arrayLB[i].y) and $(arrayLB[j])) in : \n$sols"
            end

            # non dominated solutions
            @assert !dominate(sols[i],sols[j]) "Solutions are not all nondominated, $(sols[i]) dominates $(sols[j])\n$sols"
            @assert !dominate(sols[j],sols[i]) "Solutions are not all nondominated, $(sols[j]) dominates $(sols[i])\n$sols"

            # sorted z1 increasing
            @assert sols[i].y[1] < sols[j].y[1] "The LB isn't sorted !\n$(sols)"

        end
    end
end

" `Assertion` : `correct`, `feasible`, `sorted z1 increasing`, `unicity`, `non dominated solutions` "
function DEBUG_UB(prob::BiOKP, UB::UpperBound)

    # Change type of UB.sols
    if typeof(UB.sols) != Vector{Sol}
        # convert into array (easier to work with indices)
        sols = Vector{Sol}()
        for sol in UB.sols
            push!(sols, sol)
        end
    else
        sols = UB.sols
    end
    
    # check solutions
    for i in 1:length(sols)

        # correct
        DEBUG_feasibleSolution(prob, sols[i])

        for j in (i+1):length(sols)

            # unicity
            if prob.isInteger
                @assert sols[i].y != sols[j].y "The UB contains duplicates\n$sols"
            else
                @assert differentFloatSols(sols[i],sols[j]) "The UB contains duplicates ($(arrayLB[i].y) and $(arrayLB[j])) in : \n$sols"
            end

            # non dominated solutions
            @assert !dominate(sols[i],sols[j]) "Solutions are not all nondominated, $(sols[i]) dominates $(sols[j])\n$sols"
            @assert !dominate(sols[j],sols[i]) "Solutions are not all nondominated, $(sols[j]) dominates $(sols[i])\n$sols"

            # sorted z1 increasing
            @assert sols[i].y[1] < sols[j].y[1] "The UB isn't sorted !\n$(sols)"

        end
    end

    # check segments
    if UB.explicit
        @assert (length(UB.segments.b) == 0 && length(sols) < 1) || length(UB.segments.b) == length(sols)+1 "Incorrect number of segments : $(length(UB.segments.b)) for a $(length(sols)) solutions"
        @assert (size(UB.segments.A)[1] == 0 && length(sols) < 1) || size(UB.segments.A)[1] == length(sols)+1 "Incorrect number of segments : $(size(UB.segments.A)[1]) for a $(size(UB.segments.A)[1]) solutions"

        nbConstraints = length(sols)+1
        A = zeros(Float64, nbConstraints, 2)
        b = zeros(Float64, nbConstraints)

        max1 = sols[1].y[1]
        max2 = sols[1].y[2]

        for iter in 1:length(sols)-1
            solL = sols[iter]
            solR = sols[iter+1]

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

        @assert A == UB.segments.A "Matrix A in segments is wrong : $(UB.segments.A) vs $(A)"
        @assert b == UB.segments.b "Vector b in segments is wrong : $(UB.segments.b) vs $(b)"
    end

end

" `Assertion` : `correct` "
function DEBUG_correctSolution(prob::BiOKP, sol::Sol)

    trueProfit = [0.0,0.0]
    trueWeight = 0.0
    floatNumber = false

    for i in 1:prob.nbVar
        elt = sol.x[i]
        
        if sol.isBinary
            @assert (elt != 0.0 || elt != 1.0) "Solution is not truly binary : $(sol.x), $(sol.isBinary)"
        else
            if elt != 0.0 && elt != 1.0
                floatNumber = true
            end
        end

        if sol.x[i] != 0.0
            trueWeight += sol.x[i] * round(Int8,prob.weights[i])
            trueProfit[1] += sol.x[i] * prob.profits[1,i]
            trueProfit[2] += sol.x[i] * prob.profits[2,i]
        end
    end

    sol.isBinary && @assert trueWeight == sol.w "Solution's weight is wrong : true = $trueWeight, actual = $(sol.w)"
    sol.isBinary && @assert trueProfit == sol.y "Solution's profit is wrong : true = $trueProfit, actual = $(sol.y)"
    !sol.isBinary && @assert floatNumber "Solution is TRULY binary : $floatNumber $(sol.x), $(sol.isBinary)"
end

" `Assertion` : `correct`, `feasible`, not always `binary` !"
function DEBUG_feasibleSolution(prob::BiOKP, sol::Sol)

    DEBUG_correctSolution(prob, sol)

    @assert sol.w <= prob.maxWeight "Solution is too heavy (trueWeight = $trueWeight), $(sol.w) for the problem $(prob.maxWeight)"

end

" `Assertion` : `correct`, `feasible`, `binary` "
function DEBUG_feasibleBinarySolution(prob::BiOKP, sol::Sol)
    DEBUG_feasibleSolution(prob, sol)

    @assert sol.isBinary "Solution is not binary ! $sol"
end

" `Assertion` : make sure the given `sol` isn't dominating any of the solutions in `sols` "
function DEBUG_isntDominatingAny(sol::Sol, sols::LinkedList{Sol})
    # convert into array (easier to work with indices)
    vec = Vector{Sol}()
    for sol in sols
        push!(vec, sol)
    end

    for i in 1:length(vec)
        @assert !dominate(sol, vec[i]) "$sol is dominating $(vec[i])"
    end
end

" `Assertion` : `correct`, `feasible`, `sorted increasing`, `unicity`, `non dominated solutions` "
function DEBUG_front(prob::BiOKP, front::LinkedList{Sol})
    # change linkedlist into array
    front_vec = Vector{Sol}()
    for sol in front
        push!(front_vec, sol)
    end

    DEBUG_front(prob, front_vec)
end

" `Assertion` : `correct`, `feasible`, `sorted increasing`, `unicity`, `non dominated solutions` "
function DEBUG_front(prob::BiOKP, front::Vector{Sol})
    
    for i in 1:length(front)

        DEBUG_feasibleSolution(prob, front[i])

        for j in (i+1):length(front)

            # unicity
            if prob.isInteger
                @assert front[i].y != front[j].y "The front contains duplicates\n$front"
            else
                @assert differentFloatSols(front[i],front[j]) "The front contains duplicates ($(arrayLB[i].y) and $(arrayLB[j])) in : \n$front"
            end

            # non dominated solutions
            @assert !dominate(front[i],front[j]) "Solutions are not all nondominated, $(front[i]) dominates $(front[j])\n$front"
            @assert !dominate(front[j],front[i]) "Solutions are not all nondominated, $(front[j]) dominates $(front[i])\n$front"

            # sorted
            @assert front[i].y[1] < front[j].y[1] "The lower bound isn't sorted !\n$(front)"

        end
    end
end

" `Assertion` : `correct`, `feasible`, `binary`, `left.z1 < right.z1`, `sorted increasing`, `unicity` "
function DEBUG_nadirs(prob::BiOKP,NP::Vector{PairOfSolution})

    for i in 1:(length(NP)-1)

        # feasible
        DEBUG_feasibleBinarySolution(prob, NP[i].solL)

        # pair correct : left.z1 < right.z1
        @assert NP[i].solL.y[1] < NP[i].solR.y[1] && NP[i].solL.y[2] > NP[i].solR.y[2] "Pair of nadir points is wrong\nhere : $(NP[i])\n$NP"

        # sorted and unicity
        @assert NP[i].solL.y[1] < NP[i+1].solL.y[1] "List of nadirs not sorted ! \n$NP"
    end

end

" `Assertion` : `correct`, `feasible`, `left.z1 < right.z1`, `sorted increasing`, `unicity` "
function DEBUG_pairsOfSolutions(prob::BiOKP, sols::Vector{PairOfSolution})

    for i in 1:length(sols)
        # feasible
        DEBUG_feasibleSolution(prob, sols[i].solL)
        DEBUG_feasibleSolution(prob, sols[i].solR)

        # pair correct : left.z1 < right.z1
        @assert sols[i].solL.y[1] < sols[i].solR.y[1] && sols[i].solL.y[2] > sols[i].solR.y[2] "Pair of points is wrong\nhere : $(sols[i])\n$sols"
    end

    for i in 1:length(sols)-1
        # sorted and unicity
        @assert sols[i].solL.y[1] < sols[i+1].solL.y[1] "List of points not sorted ! \n$sols"
    end
end

"""
Assertion on knownSols : `front`, `sol.x[index] == value for every sol in list`.

Assertion on previousIndicies : `sorted increasing`, `unicity`.
"""
function DEBUG_parentToChild(prob::BiOKP, PTC::ParentToChild; index = -1, value = -1)
    
    # secure knownSols
    DEBUG_front(prob, PTC.knownSols)

    if index != -1 && value != -1
        for elt in PTC.knownSols
            @assert elt.x[index] == value "The list doesn't contain the right elements"
        end
    end

    lenPrev = length(PTC.previousIndicies)

    @assert lenPrev == length(PTC.knownSols) "Known sols and previous indicies aren't of the same length"

    # secure previousIndicies
    prevArray = Vector{Int}(undef,lenPrev)
    iter = 1
    for elt in PTC.previousIndicies
        prevArray[iter] = elt
        iter += 1
    end

    for i in 1:lenPrev
        for j in (i+1):lenPrev
            @assert prevArray[i] != prevArray[j] "No unicity"
            @assert prevArray[i] < prevArray[j] "Indicies not sorted"
        end
    end
end

#---------------- PARAMETRIC RELAXATION ------------------#

" `Assertion` : `valid index`, `lambdas values sorted increasing` "
function DEBUG_lambdas(lambdas::Vector{LambdaChange})

    for i in 1:length(lambdas)

        # check if the index is a distinct tuple (not the same number) and n1 < n2
        (n1,n2) = lambdas[i].index
        @assert n1 < n2 || (n1 == 0 && n2 == 0) "The index of the lambda doesn't fit the rule n1 < n2"

        for j in (i+1):length(lambdas)

            # check if the list is sorted
            @assert lambdas[j].value <= lambdas[i].value "$i,$j,$(lambdas[j].value),$(lambdas[i].value) The lambdas are not sorted decreasing"

        end

    end

end

" `Assertion` : `nbEltNewIndicies = nbEltLambdas` "
function DEBUG_newIndices(nbEltNewIndicies::Int,nbEltLambdas::Int)
    @assert nbEltNewIndicies == nbEltLambdas "New Indicies have a different length than lambdas"
end