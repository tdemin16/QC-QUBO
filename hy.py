import dimod
import hybrid
import numpy as np

def run_annealer(theta):
    iteration = hybrid.RacingBranches(
            hybrid.InterruptableTabuSampler(),
            hybrid.EnergyImpactDecomposer(size=1)
            | hybrid.QPUSubproblemAutoEmbeddingSampler()
            | hybrid.SplatComposer()
        ) | hybrid.ArgMin()
    workflow = hybrid.LoopUntilNoImprovement(iteration, convergence=1)

    bqm = dimod.BinaryQuadraticModel({}, theta, 0, dimod.BINARY)

    init_state = hybrid.State.from_problem(bqm)
    final_state = workflow.run(init_state).result()
    response = final_state.samples.first.sample

    return np.atleast_2d(list(response.first.sample.values())).T

def to_matrix(nums):
    theta = [[0 for col in range(len(nums))] for row in range(len(nums))]
    c = sum(nums)

    for i in range(len(nums)):
        for j in range(len(nums)):
            if i != j: theta[i][j] = nums[i] * nums[j]
            else: theta[i][i] = nums[i] * (nums[i] - c)

    return np.array(theta)

def matrix_to_dict(theta):
    n = len(theta)
    d = dict()
    for i in range(n):
        for j in range(n):
            d[i, j] = theta[i][j]

    return d

def fQ(theta, sol):
    return sol.transpose() * theta * sol
            

nums = [3, 7, 4, 10, 4, 2, 10, 1, 7, 3]
theta = to_matrix(nums)

sol = run_annealer(matrix_to_dict(theta))
print(fQ(theta, sol))