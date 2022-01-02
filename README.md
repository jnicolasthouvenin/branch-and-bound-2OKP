
![baniere_b b](https://user-images.githubusercontent.com/40352310/134726925-59d748c3-0067-4d98-aecd-8c95b511a9d9.png)

# Bi-objective Branch and Bound for the monodimensional Knapsack Problem

This solver finds the set of efficient solutions for bi-objective monodimensional knapsack fully integer problems. The appellation "fully integer" denotes problems where the profits and weights are integers. The solver is based on a branch-and-bound algorithm.

The first version of this solver was implemented by Lucas Baussay, Mathilde Sedeno and myself, during a research project from february to june 2021 at the University of Nantes. Several optimizations have been added by myself between july and september 2021 to fasten the overall performances.

## Chronology

`02/2021 - 06/2021` : first version of the code (Research Project with Lucas Baussay and Mathilde Sedeno - Nantes)

`07/2021 - 09/2021` : addition of several optimizations in the computations (Research Internship at the LS2N lab - Nantes)

`12/2021 - 01/2022` : rework of the code structure and addition of a new primal heuristic for a better scaling

## How to use

1. Open Julia Repl : <code>& julia</code>

2. Import the code : <code>$ include("main.jl")</code>

3. Load an (fully integer) instance from a file in /instances : <code>$ p = FILES_readProblem(fileName)</code>, or a randomly generated (fully integer) instance : <code>$ p = BiOKP(nbVars)</code>

4. Solve it using the command : <code>$ main(p)</code>

5. Print the computation times of each part of the code : <code>$ show(to)</code>

## Doc

<table>
<thead>
<tr>
<th>Code file</th>
<th>Label</th>
</tr>
</thead>
<tbody>
<tr>
<td>GeneticAlgorithms</td>
<td>Module behaving as primal heuristics</td>
</tr>
<tr>
<td>GAInterface.jl</td>
<td>Interface for communicating with the module GeneticAlgorithms</td>
</tr>
<tr>
<td>assignment.jl</td>
<td>Functions for handling the modification of the assignements in the branch and bound</td>
</tr>
<tr>
<td>common.jl</td>
<td>Common operations made in the code</td>
</tr>
<tr>
<td>debug.jl</td>
<td>Assertions on the variables used in the program. Helps locating errors right away</td>
</tr>
<tr>
<td>files.jl</td>
<td>Functions that write and read 2OKP instances on files</td>
</tr>
<tr>
<td>fractionsStruct.jl</td>
<td>Like rational type but using float numbers for numerator and denumerator (specific to Parametric Relaxation)</td>
</tr>
<tr>
<td>lowerBound.jl</td>
<td>Working with the LowerBound</td>
</tr>
<tr>
<td>main.jl</td>
<td>Main file</td>
</tr>
<tr>
<td>parametricLinearRelax.jl</td>
<td>Implementation of the parametric method for linear relaxation</td>
</tr>
<tr>
<td>solve1OKP.jl</td>
<td>Functions that solve 1OKP problems</td>
</tr>
<tr>
<td>structs.jl</td>
<td>Declaration of the original types used in the algorithm</td>
</tr>
<tr>
<td>subProblemStudy.jl</td>
<td>Functions indicating if a subproblem can be fathomed and giving additional information</td>
</tr>
<tr>
<td>tests.jl</td>
<td>Tests on the code</td>
</tr>
<tr>
<td>upperBound.jl</td>
<td>Upper bound computation</td>
</tr>
<tr>
<td>utils.jl</td>
<td>Basic utilitary functions</td>
</tr>

</tbody>
</table>
