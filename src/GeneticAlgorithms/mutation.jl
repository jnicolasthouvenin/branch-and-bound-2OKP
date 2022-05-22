
function MUTATION_flip!(offspring::Sol)

    debug && DEBUG_correct_solution(offspring)

    mut_ind = rand(1:offspring.prob.n) # mutation index
    offspring.x[mut_ind] == 0 ? SOLUTION_add_item!(offspring, mut_ind) : SOLUTION_remove_item!(offspring, mut_ind) # flip

    debug && DEBUG_correct_solution(offspring)

end
