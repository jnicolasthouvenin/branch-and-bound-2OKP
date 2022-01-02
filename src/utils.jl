
""" Basic utilitary functions """

# All functions are public

" True if `sol1` and `sol2` are different solutions, false otherwise. "
function differentFloatSols(sol1::Sol,sol2::Sol)
    return abs(sol1.y[1]-sol2.y[1]) > 0.000001 || abs(sol1.y[2]-sol2.y[2]) > 0.000001
end

# NOTE : for domination tests, one could overwrite the operator >= of julia
# and test dominance using sol1 >= sol2

" True if the given solution `sol1` dominates the solution `sol2`. "
function dominate(sol1::Sol, sol2::Sol)
	return (sol1.y[1] >= sol2.y[1] && sol1.y[2] > sol2.y[2]) || (sol1.y[1] > sol2.y[1] && sol1.y[2] >= sol2.y[2])
end

" True if the given solution `sol` dominates the nadirPoint `nadirPoint`. "
function dominate(sol::Sol, nadirPoint::PairOfSolution)
	return sol.y[1] >= (nadirPoint.solL.y[1] + 1) && sol.y[2] >= (nadirPoint.solR.y[2]+1)
end
