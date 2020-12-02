import dimod
import hybrid
import os
import sys
import dwave_networkx as dnx
import networkx as nx
import scipy


def run_annealer(theta, iteration, workflow):
    # Build the QUBO problem
    bqm = dimod.BinaryQuadraticModel({}, theta, 0, dimod.SPIN)

    # Solve
    init_state = hybrid.State.from_problem(bqm)
    final_state = workflow.run(init_state).result()
    solution = final_state.samples.first.sample

    return solution

def chimera(r,c):
    G = dnx.chimera_graph(r, c)
    tmp = nx.to_dict_of_lists(G)
    n = len(tmp)
    rows = []
    cols = []
    for i in range(n):
        rows.append(i)
        cols.append(i)
        for element in tmp[i]:
            rows.append(i)
            cols.append(element)
    
    return list(zip(rows, cols))


def main():
    n = sys.stdin.read(6)
    n = int(n.split('\x00', 1)[0])
    
    simulation = sys.stdin.read(2)
    simulation = int(simulation.split('\x00', 1)[0])

    row = 16 * 8
    if(n >= row): n_cols = 16
    else: n_cols = n / 8
    n_rows = int(n / row)
    if n % row != 0:
        n_rows += 1

    A = chimera(n_rows, n_cols)

    for r, c in A:
        msg = str(r)
        msg = ('0' * (4 - len(msg)) + msg).encode()
        os.write(1, msg)

        msg = str(c)
        msg = ('0' * (4 - len(msg)) + msg).encode()
        os.write(1, msg)
        pass
    msg = ("####").encode()
    os.write(1, msg)
    os.write(1, msg)


    if(simulation == 0):
        theta = dict()
        i = 0
        iteration = hybrid.RacingBranches(
            hybrid.InterruptableTabuSampler(),
            hybrid.EnergyImpactDecomposer(size=1)
            | hybrid.QPUSubproblemAutoEmbeddingSampler()
            | hybrid.SplatComposer()
        ) | hybrid.ArgMin()

        workflow = hybrid.LoopUntilNoImprovement(iteration, convergence=1)

        end = False
        while end == False:
            x = sys.stdin.read(100)
            if(x[0] != "#" and x[0:3] != "END"):  # retrieving problem from pipe
                x = x.split('\x00', 1)[0]
                if i == 0:
                    r = int(x)
                elif i == 1:
                    c = int(x)
                else:
                    theta[r, c] = float(x)
                    pass

                i = (i + 1) % 3

            elif(x[0] == "#"):  # compute problem
                i = 0
                l = run_annealer(theta, iteration, workflow)

                for j in range(len(l)):
                    if(l[j] == 1):
                        msg = ("+" + str(l[j])).encode()
                    else:
                        l[j] = 1
                        msg = ("-" + str(l[j])).encode()
                    os.write(1, msg)  # write solution on pipe
                    pass

                theta.clear()  # clear dictionary
                pass
            
            else:
                end = True
                pass
            pass


if __name__ == '__main__':
    main()
