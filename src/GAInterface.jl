
""" Provide functions to communicate with the module `GeneticAlgorithms` """

# Public functions :
# - GAINTERFACE_exportProb : Convert `BiOKP` to `_MOMKP`
# - GAINTERFACE_exportLB : Convert `LowerBound` into a `front` for `GeneticAlgorithms`
# - GAINTERFACE_importLB : Convert a `front` for `GeneticAlgorithms` to a `LowerBound`

# Export functions convert types of the code to the types of the module GeneticAlgorithms
# Import functions convert types of the module GeneticAlgorithms to the types of the code

#-------------------------------------------------------------------#
#------------------------ PUBLIC FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

"""
Convert `BiOKP` to `_MOMKP`
"""
function GAINTERFACE_exportProb(prob::BiOKP)
    @assert prob.isInteger "A _MOMKP object can only be integer"

    function weights_int(weights::Vector{Float64})
        W = Matrix{Int}(undef, 1, length(weights))
        for i in 1:length(weights)
            W[1,i] = Int(weights[i])
        end
        return W
    end

    return GeneticAlgorithms._MOMKP("myMOMKP",
                prob.nbVar,
                1,
                round.(Int, prob.profits),
                weights_int(prob.weights),
                [Int(prob.maxWeight)]
            )
end

"""
Convert `LowerBound` into a `front` for `GeneticAlgorithms`
"""
function GAINTERFACE_exportLB(prob::BiOKP, LB::LowerBound)
    n = length(LB.sols)
    LB_vec = Vector{GeneticAlgorithms.Sol}(undef, n)
    
    LB_aux = LB.sols # (LinkedList{Sol})
    iter = 1
    while LB_aux.tail != nil(Sol)
        # fill the LB (array)
        LB_vec[iter] = GAINTERFACE_exportSol(prob, LB_aux.head)
        
        LB_aux = LB_aux.tail # go next
        iter += 1
    end
    LB_vec[iter] = GAINTERFACE_exportSol(prob, LB_aux.head)

    return LB_vec
end

"""
Convert a `front` for `GeneticAlgorithms` to a `LowerBound`
"""
function GAINTERFACE_importLB(prob::BiOKP, front::Vector{GeneticAlgorithms.Sol})
    sort!(front, by = sol -> sol.z[1], rev = true) # the construction of the LB will reverse it again

    LB = LowerBound()

    LB.sols = cons(GAINTERFACE_importSol(front[1]), LB.sols)

    for iter in 2:length(front)
        LB.sols = cons(GAINTERFACE_importSol(front[iter]), LB.sols)
    end

    CONFIG.debug && DEBUG_LB(prob, LB)

    return LB
end

#-------------------------------------------------------------------#
#----------------------- PRIVATE FUNCTIONS -------------------------#
#-------------------------------------------------------------------#

"""
Convert type `Sol` to type `GeneticAlgorithms.Sol`
"""
function GAINTERFACE_exportSol(prob::BiOKP, sol::Sol)
    @assert sol.isBinary "GeneticAlgorithms module doesn't support non binary solutions"
    return GeneticAlgorithms.Sol(GAINTERFACE_exportProb(prob), round.(Int, sol.x), round.(Int, sol.y), [round(Int, sol.w)])
end

"""
Convert type `GeneticAlgorithms.Sol` to type `Sol`
"""
function GAINTERFACE_importSol(sol::GeneticAlgorithms.Sol)
    return Sol(sol.x, sol.z, sol.w[1], true)
end