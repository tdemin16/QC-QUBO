import dimod
import hybrid
import time
import numpy as np
import os
import sys


def dict_to_vector(dic):
    n = len(dic)
    vector = []
    for i in range(n):
        vector.append(dic[i])

    return vector


def run_annealer(Q):
    # Build the QUBO problem
    bqm = dimod.BinaryQuadraticModel({}, Q, 0, dimod.SPIN)

    start = time.time_ns()
    # Define the workflow
    iteration = hybrid.RacingBranches(
        hybrid.InterruptableTabuSampler(),
        hybrid.EnergyImpactDecomposer(size=1)
        | hybrid.QPUSubproblemAutoEmbeddingSampler()
        | hybrid.SplatComposer()
    ) | hybrid.ArgMin()
    end = time.time_ns() - start
    f = open("out.txt", "a")
    f.write(str(end))
    f.close()

    workflow = hybrid.LoopUntilNoImprovement(iteration, convergence=1)

    # Solve
    init_state = hybrid.State.from_problem(bqm)
    final_state = workflow.run(init_state).result()
    solution = final_state.samples.first.sample

    return dict_to_vector(solution)


def main():
    Q = dict()
    i = 0
    finish = False

    while True:
        x = sys.stdin.read(100)
        if(x[0] != "#"):
            x = x.split('\x00',1)[0]
            if i == 0:
                r = int(x)
            elif i == 1:
                c = int(x)
            else:
                Q[r, c] = float(x)
                pass

            i = (i + 1) % 3

        else:
            i = 0
            l = run_annealer(Q)
            for j in l:
                if(j == 1):
                    msg = ("+" + str(j)).encode()
                else:
                    j = 1
                    msg = ("-" + str(j)).encode()
                os.write(1, msg)
                pass

            Q.clear()
            pass

        pass


if __name__ == '__main__':
    main()
