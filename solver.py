import dimod
import hybrid
import os
import sys
import dwave_networkx as dnx
import networkx as nx


def run_annealer(theta, iteration, workflow):
    # Build the QUBO problem
    bqm = dimod.BinaryQuadraticModel({}, theta, 0, dimod.SPIN)

    # Solve
    init_state = hybrid.State.from_problem(bqm)
    final_state = workflow.run(init_state).result()
    solution = final_state.samples.first.sample

    return solution


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

    chimera_topology = dnx.chimera_graph(n_rows, n_cols)
    A = nx.adjacency_matrix(chimera_topology)
    r, c = A.nonzero()
    for r, c in zip(r, c):
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
            x = sys.stdin.read(4096)
            if(x[0] != "#"):  # retrieving problem from pipe
                x = x.split(',')
                for k in x:
                    k = k.split('\x00', 1)[0]
                    if(i == 0):
                        r = int(k)
                    elif(i == 1):
                        c = int(k)
                    else: 
                        theta[r, c] = float(k)

                    i = (i + 1) % 3
                
            elif(x[0] == "#"):
                i = 0
                l = run_annealer(theta, iteration, workflow)
                #l = {0:-1, 1:1, 2:-1, 3:1, 4:1, 5:-1, 6:1, 7:-1, 8:1, 9:1, 10:-1, 11:1, 12:-1, 13:1, 14:1, 15:-1}
                theta.clear()
                msg = ""
                size = 0
                for j in l:
                    msg = msg + str(l[j]) + ","
                    size += len(str(l[j])) + 1

                    if(size > 4000):
                        msg = (msg + ('\0' * (4096 - size))).encode()
                        os.write(1, msg)
                        msg = ""
                        size = 0
                        pass
                    pass

                msg = (msg + ('\0' * (4096 - size))).encode()
                os.write(1, msg)
                pass 
            else:
                end = True   
            pass


if __name__ == '__main__':
    main()
