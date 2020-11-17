import dimod
import hybrid
import time
import numpy as np
import sys

def dict_to_vector(dic):
    n = len(dic)
    vector = []
    for i in range(n):
        vector.append(dic[i])

    return vector

def run_annealer(Q, ret_dict = False):    
    # Build the QUBO problem
    if isinstance(Q, dict):
        bqm = dimod.BinaryQuadraticModel({}, Q, 0, dimod.SPIN)
    elif isinstance(Q, np.ndarray):
        new_Q = matrix_to_dict(Q)
        bqm = dimod.BinaryQuadraticModel({}, new_Q, 0, dimod.SPIN)
    else:
        print("Q type incorrect")
        raise TypeError
    
    # Define the workflow
    iteration = hybrid.RacingBranches(
        hybrid.InterruptableTabuSampler(),
        hybrid.EnergyImpactDecomposer(size=1)
        | hybrid.QPUSubproblemAutoEmbeddingSampler()
        | hybrid.SplatComposer()
    ) | hybrid.ArgMin()
    
    workflow = hybrid.LoopUntilNoImprovement(iteration, convergence=1)

    # Solve
    init_state = hybrid.State.from_problem(bqm)
    final_state = workflow.run(init_state).result()
    solution = final_state.samples.first.sample
    
    return dict_to_vector(solution)

def main(n):
    j_max = 0
    j = 0
    Q = dict()
    for i in range(n):
        j_max += 1
        while j < j_max:
            Q[i,j] = np.random.randint(low=-10, high=10)
            Q[j,i] = Q[i,j]
            j += 1
        j = 0
    l = run_annealer(Q)
    print(l)

if __name__ == '__main__':
    main(8)