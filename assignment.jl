
""" Functions for handling the modification of the assignements in the branch and bound """

# All functions are public

"""
Assign the next variable to one.
"""
function ASSIGNMENT_addVar!(assignment::Assignment, prob::BiOKP)
    
    @timeit to "ASSIGNMENT_addVar!" begin

	assignment.lastAssigned += 1
	assignment.assign[assignment.lastAssigned] = 1

	assignment.weight += prob.weights[assignment.lastAssigned]
	
	assignment.profit += prob.profits[1:end, assignment.lastAssigned]

    end # TimerOutputs
end

"""
Assign the current variable to zero.

The parameter `canAddVar` indicates if the variable has been assigned to one in the previous iteration. If `canAddVar` is true, we have to remove the weight and profit of the variable from the assignment.
"""
function ASSIGNMENT_removeVar!(assignment::Assignment, prob::BiOKP, canAddVar::Bool)
    
    @timeit to "ASSIGNMENT_removeVar!" begin
	
    if !canAddVar
		assignment.lastAssigned += 1
	else
		assignment.weight -= prob.weights[assignment.lastAssigned]
		assignment.profit -= prob.profits[1:end, assignment.lastAssigned]
	end

	assignment.assign[assignment.lastAssigned] = 0

    end # TimerOutputs
end

"""
Update the `assignment` to match the version of the parent (unassigned the current variable).

NOTE : no need to take of the weight of object ... because the backtrack always happens after ASSIGNMENT_removeVar!
"""
function ASSIGNMENT_backtrack!(assignment::Assignment)
    
    @timeit to "ASSIGNMENT_backtrack!" begin

	assignment.assign[assignment.lastAssigned] = -1
	assignment.lastAssigned -= 1

    end # TimerOutputs
end