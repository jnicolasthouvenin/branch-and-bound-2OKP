
""" Functions that solve 1OKP problems """

# All functions are public

" struct to give to the Combo C library "
struct Combo_item
    p::Clonglong
    w::Clonglong
    x::Cint
    i::Cint
end

"""
Returns the optimal solution for the scalarized exact problem computed by the Combo solver.

`λ` is the scalarization vector.

`assignment` is the current assignment to take into account.
"""
function SOLVE1OKP_combo(prob::BiOKP, λ::Vector{Float64}, assignment::Assignment)

    # Handling extreme senarios
	if prob.nbVar == assignment.lastAssigned # no free variables
        returnedSol = Sol(Float64.(assignment.assign), assignment.profit, assignment.weight, true)

        CONFIG.debug && DEBUG_feasibleBinarySolution(prob,returnedSol)

		return returnedSol
	else
        # check if we are in the case where all the objects can be added or if no object can be added at all

		objectsCanBeAdded = false
        iter = assignment.lastAssigned+1
		while !objectsCanBeAdded && iter <= prob.nbVar
			objectsCanBeAdded = objectsCanBeAdded || prob.weights[iter] <= (prob.maxWeight - assignment.weight)
            iter += 1
		end

		if !objectsCanBeAdded # we can't add anymore objects

            sol = Sol(append!(Float64.(assignment.assign[1:assignment.lastAssigned]), zeros(Float64, prob.nbVar-assignment.lastAssigned)), assignment.profit, assignment.weight, true)

            CONFIG.debug && DEBUG_feasibleBinarySolution(prob,sol)

			return sol

		elseif sum(prob.weights[(assignment.lastAssigned+1):end]) <= (prob.maxWeight - assignment.weight) # we can add all the objects

            sol = Sol(append!(Float64.(assignment.assign[1:assignment.lastAssigned]), ones(Float64, prob.nbVar-assignment.lastAssigned)), assignment.profit + [sum(prob.profits[iter, (assignment.lastAssigned+1):end]) for iter=1:prob.nbObj], assignment.weight + sum(prob.weights[(assignment.lastAssigned+1):end]), true)

            CONFIG.debug && DEBUG_feasibleBinarySolution(prob,sol)

			return sol

		end
	end

    @timeit to "SOLVE1OKP_combo" begin

    weightedProfits = Vector{Float64}(undef,prob.nbVar)
    for var in 1:prob.nbVar
        weightedProfits[var] = (λ[1]*prob.profits[1,var] + λ[2]*prob.profits[2,var])
    end

	items = [Combo_item(weightedProfits[iter], prob.weights[iter], 0, iter) for iter = (assignment.lastAssigned+1):prob.nbVar]

	z = ccall((:solve, comboPath),Clonglong,(Ref{Combo_item}, Cint, Clonglong, Clonglong, Clonglong),items, (prob.nbVar - assignment.lastAssigned), (prob.maxWeight - assignment.weight), 0, 0)

    if z == 0
        z = dot(weightedProfits[(assignment.lastAssigned+1):end], broadcast(it->it.x, items))
    end

    x = Float64.(append!(assignment.assign[1:assignment.lastAssigned], falses(prob.nbVar-assignment.lastAssigned)))
    for it in items
    	x[it.i] = Float64(it.x)
    end

    returnedSol = Sol(x, assignment.profit + [sum(x[(assignment.lastAssigned+1):end] .* prob.profits[iter, (assignment.lastAssigned+1):end]) for iter = 1:prob.nbObj],assignment.weight + sum(prob.weights[(assignment.lastAssigned+1):end] .* x[(assignment.lastAssigned+1):end]), true)

    end # TimerOutput

    CONFIG.debug && DEBUG_feasibleBinarySolution(prob,returnedSol)

    return returnedSol
end

"""
Returns the optimal solution for the scalarized and relaxed problem (continuous relaxation).

`λ` is the scalarization vector.

`assignment` is the current assignment to take into account.
"""
function SOLVE1OKP_linear(prob::BiOKP, λ::Vector{Float64}, assignment::Assignment)

    @timeit to "SOLVE1OKP_linear" begin

    utilities = Vector{Float64}(undef,prob.nbVar-assignment.lastAssigned)

    for i in 1:prob.nbVar-assignment.lastAssigned
        utilities[i] = (λ[1] * prob.profits[1,i+assignment.lastAssigned] + λ[2] * prob.profits[2,i+assignment.lastAssigned]) ./ prob.weights[i+assignment.lastAssigned]
    end

    permList = sortperm(utilities, rev = true)

	nbUnassignedVars = prob.nbVar-assignment.lastAssigned

	weight = assignment.weight
	profit = assignment.profit[1:end]

    x = Float64.(assignment.assign[1:end])

    for iter = 1:nbUnassignedVars
        permList[iter] += assignment.lastAssigned
        x[assignment.lastAssigned + iter] = 0.
    end

    iter = 1
    isBinary = true
    while iter <= nbUnassignedVars && weight + prob.weights[permList[iter]] <= prob.maxWeight
        weight += prob.weights[permList[iter]]
        x[permList[iter]] = 1.
        profit[1] += prob.profits[1,permList[iter]]
        profit[2] += prob.profits[2,permList[iter]]

        iter += 1
    end

    if iter <= nbUnassignedVars && (prob.maxWeight - weight) > 0.
        x[permList[iter]] = (prob.maxWeight - weight) / prob.weights[permList[iter]]
        profit += x[permList[iter]] * prob.profits[1:end, permList[iter]]
        weight += x[permList[iter]] * prob.weights[permList[iter]]
        isBinary = false
    end

    sol = Sol(x, profit, weight, isBinary)

    end # TimerOutput

    CONFIG.debug && DEBUG_feasibleSolution(prob,sol)

	return sol
end

"""
Returns the optimal solution for the scalarized exact problem computed by the solver GLPK (using JuMP).

`λ` is the scalarization vector.

`assignment` is the current assignment to take into account.
"""
function SOLVE1OKP_GLPK(prob::BiOKP, λ::Vector{Float64}, assignment::Assignment)

    @timeit to "SOLVE1OKP_GLPK" begin

	model = Model(GLPK.Optimizer)
	x = @variable(model, x[1:(prob.nbVar-assignment.lastAssigned)], Bin)
	@constraint(model, Weights, sum(x .* prob.weights[(assignment.lastAssigned+1):end]) + assignment.weight <= prob.maxWeight)
	@objective(model, Max, λ[1] * sum(x .* prob.profits[1, (assignment.lastAssigned+1):end]) + λ[2] * sum(x .* prob.profits[2, (assignment.lastAssigned+1):end]))

	optimize!(model)

	X = append!(assignment.assign[1:assignment.lastAssigned],(Float64).(value.(x)))

	termStatus = termination_status(model)

    end # TimerOutput

	if termStatus == MOI.OPTIMAL
        sol = Sol(X, [sum(X .* prob.profits[1,1:end]), sum(X .* prob.profits[2,1:end])], sum(X .* prob.weights), true)

        CONFIG.debug && DEBUG_feasibleBinarySolution(prob,sol)

        return sol

	elseif termStatus == MOI.INFEASIBLE
        error("GLPK is stating that the problem is infeasible")
	else
		error("MOI or GLPK aren't working as they should be")
	end
end