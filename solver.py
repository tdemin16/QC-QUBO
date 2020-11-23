import dimod
import hybrid
import os
import sys


def dict_to_vector(dic):
    n = len(dic)
    vector = []
    for i in range(n):
        vector.append(dic[i])

    return vector


def run_annealer(Q, iteration, workflow):
    # Build the QUBO problem
    bqm = dimod.BinaryQuadraticModel({}, Q, 0, dimod.SPIN)

    # Solve
    init_state = hybrid.State.from_problem(bqm)
    final_state = workflow.run(init_state).result()
    solution = final_state.samples.first.sample

    return solution


def main():
    Q = dict()
    i = 0
    finish = False

    iteration = hybrid.RacingBranches(
        hybrid.InterruptableTabuSampler(),
        hybrid.EnergyImpactDecomposer(size=1)
        | hybrid.QPUSubproblemAutoEmbeddingSampler()
        | hybrid.SplatComposer()
    ) | hybrid.ArgMin()

    workflow = hybrid.LoopUntilNoImprovement(iteration, convergence=1)

    while True:
        x = sys.stdin.read(100)
        if(x[0] != "#"): # retrieving problem from pipe
            x = x.split('\x00',1)[0]
            if i == 0:
                r = int(x)
            elif i == 1:
                c = int(x)
            else:
                Q[r, c] = float(x)
                pass

            i = (i + 1) % 3

        else: #compute problem
            i = 0
            l = run_annealer(Q, iteration, workflow)
            
            for j in range(len(l)):
                if(l[j] == 1):
                    msg = ("+" + str(l[j])).encode()
                else:
                    l[j] = 1
                    msg = ("-" + str(l[j])).encode()
                os.write(1, msg) # write solution on pipe
                pass

            Q.clear() # clear dictionary
            pass

        pass


if __name__ == '__main__':
    main()
