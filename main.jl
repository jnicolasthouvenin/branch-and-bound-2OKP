
""" Main file """

# All functions are public

#--------------------------- LOAD LIBRARIES ---------------------------#

println("\n\$ Loading libraries")

using Printf
@printf("%-6s %-7s %-21s %s\n", "[INFO]", "using", "Printf", "")
@printf("%-6s %-7s %-21s %s\n", "[INFO]", "using", "TimerOutputs", "(measure times)")
using TimerOutputs
@printf("%-6s %-7s %-21s %s\n", "[INFO]", "using", "DataStructures", "(LinkedList)")
using DataStructures
@printf("%-6s %-7s %-21s %s\n", "[INFO]", "using", "Random", "")
using Random
@printf("%-6s %-7s %-21s %s\n", "[INFO]", "using", "JuMP", "(FrameWork for GLKP solver)")
using JuMP
@printf("%-6s %-7s %-21s %s\n", "[INFO]", "using", "GLPK", "(Solver)")
using GLPK
@printf("%-6s %-7s %-21s %s\n", "[INFO]", "using", "GeneticAlgorithms", "(Primal Heuristic)")
include("GeneticAlgorithms/main.jl")

#--------------------------- GLOBAL VARIABLES ---------------------------#

const comboPath = joinpath(@__DIR__,"..","deps","libcombo.so")

@enum FirstDicho JUMP COMBO
@enum MethodUB EXACT_JUMP EXACT_COMBO RELAX_LIN_CLASSIC RELAX_LIN_SPEED_UP

# Algorithm components selected for the execution
mutable struct Components
    firstDicho::FirstDicho
    methodUB::MethodUB
    nadirsShift::Bool
    primal_heuristic::Bool

    Components() = new(COMBO, EXACT_COMBO, true, true) # best config possible (see experimentations)
end

# Configuration
mutable struct Config
    debug::Bool
    # verbose
    # ...

    Config() = new(true) # by default we are in debug mode
end

CONFIG     = Config()
COMPONENTS = Components()
to         = TimerOutput() # stores the time spent on the different functions

#--------------------------- CODE ---------------------------#

println("\n\$ Loading source files")

# Elementary
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "fractionsStruct", "(parametricRelax)")
include("fractionsStruct.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "structs", "(core)")
include("structs.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "GAInterface", "(GeneticAlgorithms)")
include("GAInterface.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "utils", "(core)")
include("utils.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "assignment", "(core)")
include("assignment.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "debug", "(core)")
include("debug.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "files", "(core)")
include("files.jl")

# Algorithms
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "solve1OKP", "(core)")
include("solve1OKP.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "common", "(core)")
include("common.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "lowerBound", "(core)")
include("lowerBound.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "upperBound", "(core)")
include("upperBound.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "subProblemStudy", "(core)")
include("subProblemStudy.jl")
@printf("%-6s %-7s %-21s %-s\n", "[INFO]", "include", "parametricLinearRelax", "(parametricRelax)")
include("parametricLinearRelax.jl")

println("")

"""
Returns the new nadir points to study and the time spent on computation. Modifies the LB in place.
"""
function branchAndBound!(prob::BiOKP, LB::LowerBound, assignment::Assignment, nadirsToStudy::Vector{PairOfSolution}; verbose = false, LSU::LinearSpeedUp = LinearSpeedUp(), parentToChild = ParentToChild(), num::Int = 1, iterator = nothing, timeMax = nothing, start = nothing)

    CONFIG.debug && DEBUG_LB(prob, LB)
    CONFIG.debug && DEBUG_nadirs(prob, nadirsToStudy)
    CONFIG.debug && DEBUG_parentToChild(prob, parentToChild)

    # extreme scenario : exceed time limit
	if timeMax != nothing && time() - start > timeMax
		return nadirsToStudy, false
	end

    @timeit to "b&b" begin

	testTime = false

	iterator != nothing && (iterator.value += 1)

    # Study the subproblem
    if COMPONENTS.methodUB == RELAX_LIN_SPEED_UP
        prunedType, dominatedNadirs, nonDomNadirs, UB = SPS_parametric(prob, assignment, nadirsToStudy, LSU, parentToChild = parentToChild)
    else
        prunedType, dominatedNadirs, nonDomNadirs, UB = SPS_partially_implicit(prob, assignment, nadirsToStudy, parentToChild = parentToChild)
    end
    
    #prunedType, dominatedNadirs, nonDomNadirs, UB = SPS_subProblemStudy(prob, assignment, nadirsToStudy, parentToChild = parentToChild, LSU = LSU)

    # update the LB if necessary
    if prunedType == NONE || prunedType == OPTIMALITY
		newNadirsToStudy = LB_update!(prob, LB, dominatedNadirs, UB)
    end

    end # TimerOutput

	if prunedType == NONE && assignment.lastAssigned < prob.nbVar
		canAddVar = (assignment.weight + prob.weights[assignment.lastAssigned + 1] <= prob.maxWeight)

        # if the variable can be added without breaking the capacity constraint
		if canAddVar
			ASSIGNMENT_addVar!(assignment, prob)
			parentToChild = UB_computeParentToChild(prob, UB, assignment.lastAssigned, 1)

            newNadirsToStudy, testTime = branchAndBound!(prob,
                                                LB,												
												assignment,
												newNadirsToStudy,
                                                LSU = LSU,
												parentToChild = parentToChild,
												num = num + 1,
												iterator = iterator,
												timeMax = timeMax,
												start = start)
		end

		ASSIGNMENT_removeVar!(assignment, prob, canAddVar)
		parentToChild = UB_computeParentToChild(prob, UB, assignment.lastAssigned, 0)

        newNadirsToStudy, testTime = branchAndBound!(prob,
											LB,
											assignment,
											newNadirsToStudy,
                                            LSU = LSU,
											parentToChild = parentToChild,
											num = num + 1,
											iterator = iterator,
											timeMax = timeMax,
											start = start)

		ASSIGNMENT_backtrack!(assignment)

        CONFIG.debug && DEBUG_LB(prob, LB)
        CONFIG.debug && DEBUG_nadirs(prob, newNadirsToStudy)
        CONFIG.debug && DEBUG_nadirs(prob, nonDomNadirs)
        CONFIG.debug && DEBUG_parentToChild(prob, parentToChild)

		return mergeNadirs(prob, newNadirsToStudy, nonDomNadirs), (timeMax == nothing || time() - start < timeMax)
	elseif prunedType == OPTIMALITY

        CONFIG.debug && DEBUG_LB(prob, LB)
        CONFIG.debug && DEBUG_nadirs(prob, newNadirsToStudy)
        CONFIG.debug && DEBUG_nadirs(prob, nonDomNadirs)
        CONFIG.debug && DEBUG_parentToChild(prob, parentToChild)

		return mergeNadirs(prob,newNadirsToStudy, nonDomNadirs), (timeMax == nothing || time() - start < timeMax)
	else

        CONFIG.debug && DEBUG_LB(prob, LB)
        CONFIG.debug && DEBUG_nadirs(prob, nadirsToStudy)
        CONFIG.debug && DEBUG_parentToChild(prob, parentToChild)

		return nadirsToStudy, (timeMax == nothing || time() - start < timeMax)
	end
end

"""
Returns the `LB` for the given problem. Also specifies the number of subproblems studied during execution.

The problem `prob` is assumed to have integer constants and exactly two objectives.

`WARNING` : solving medium to large instances can lead to excessively long resolutions.

# EXAMPLES

```jldoctest
julia> p = FILES_readProblem("Inst_10_1")
BiOKP(10, 2, [50.0 50.0 … 22.0 7.0; 42.0 20.0 … 1.0 16.0], [13.0, 41.0, 29.0, 27.0, 21.0, 30.0, 43.0, 42.0, 14.0, 40.0], 120.0, true)

julia> main(p)
(31, 
[93.0, 157.0]
[136.0, 142.0]
[148.0, 141.0]
[153.0, 137.0]
[164.0, 123.0]
[177.0, 118.0]
[189.0, 117.0]
[191.0, 87.0]
)
```

# DIDACTIC

```jldoctest
julia> main()
(19, 
[10.0, 15.0]
[13.0, 14.0]
[21.0, 12.0]
[22.0, 10.0]
[30.0, 7.0]
)
```
"""
function main(prob::BiOKP; timeMax = nothing, start = nothing)
    
    reset_timer!(to) # TimerOutput

	@assert prob.nbObj == 2 "This Branch and Bound supports only a Bio-objective problem"

    # secure resolution parameters
    if !prob.isInteger
        @assert false "The current algorithm only support instances with integer parameters"
        if COMPONENTS.methodUB == EXACT_COMBO
            COMPONENTS.methodUB = EXACT_JUMP
        end
        COMPONENTS.nadirsShift = false
    end

	assignment = Assignment(prob)

	compt = Iterator()

    LB = LB_dichoSearch(prob, assignment)

    if COMPONENTS.primal_heuristic && prob.nbVar >= 50 # don't call the GAs on very small instances
        momkp        = GAINTERFACE_exportProb(prob)
        seeding_sols = GAINTERFACE_exportLB(prob, LB)
        # by default time spent on GAs is equal to prob.nbVar/10
        front = GeneticAlgorithms.solve(momkp, max_t = Float64(prob.nbVar)/10, seeding_sols = seeding_sols)
        LB           = GAINTERFACE_importLB(prob, front)
    end

	nadirs = LB_getNadirs(prob, LB)

    if COMPONENTS.methodUB == RELAX_LIN_SPEED_UP
        LSU = initialiseLinearSpeedUp(prob)
    else
        LSU = LinearSpeedUp()
    end

	branchAndBound!(prob, LB, assignment, nadirs, LSU = LSU, parentToChild = ParentToChild(), iterator = compt, timeMax = timeMax, start = start)

    CONFIG.debug && DEBUG_LB(prob,LB)

	return compt.value, LB
end

function main(fname::String)
	prob = BiOKP(fname)
	return main(prob)
end

function main()
	prob = BiOKP()
	return main(prob)
end
