
function CROSSOVER_binary_uniform(p1::Sol, p2::Sol)

    debug && DEBUG_feasible_solution(p1)
    debug && DEBUG_feasible_solution(p2)

    momkp = p1.prob
    n = momkp.n

    uniform_sample = rand(1:2, n)

    offspring = SOLUTION_empty_solution(momkp)

    for idx_var in 1:n
        # get the value of the gene
        value = p1.x[idx_var]
        if uniform_sample[idx_var] == 2
            value = p2.x[idx_var]
        end

        # if the item has been added, modify the solution
        if value == 1
            SOLUTION_add_item!(offspring, idx_var)
        end
    end

    debug && DEBUG_correct_solution(offspring)

    return offspring

end

function CROSSOVER_one_point(p1::Sol, p2::Sol)

    debug && DEBUG_feasible_solution(p1)
    debug && DEBUG_feasible_solution(p2)

    momkp = p1.prob
    n = momkp.n

    point = rand(1:n)

    offspring = SOLUTION_empty_solution(momkp)

    for idx_var in 1:n
        # get the value of the gene
        if idx_var <= point
            value = p1.x[idx_var]
        else
            value = p2.x[idx_var]
        end

        # if the item should be been added, modify the solution
        if value == 1
            SOLUTION_add_item!(offspring, idx_var)
        end
    end

    debug && DEBUG_correct_solution(offspring)

    return offspring

end

function CROSSOVER_two_points(p1::Sol, p2::Sol)

    debug && DEBUG_feasible_solution(p1)
    debug && DEBUG_feasible_solution(p2)

    momkp = p1.prob
    n = momkp.n

    point_1 = rand(1:n-1)
    point_2 = rand(point_1+1:n)

    offspring = SOLUTION_empty_solution(momkp)

    for idx_var in 1:n
        # get the value of the gene
        if idx_var <= point_1
            value = p1.x[idx_var]
        elseif idx_var <= point_2
            value = p2.x[idx_var]
        else
            value = p1.x[idx_var]
        end

        # if the item should be been added, modify the solution
        if value == 1
            SOLUTION_add_item!(offspring, idx_var)
        end
    end

    debug && DEBUG_correct_solution(offspring)

    return offspring

end