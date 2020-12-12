import dimod
import os
import sys
import dwave_networkx as dnx
import networkx as nx
import scipy
import signal
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite


def handler(signum, frame):
    exit(0)

def run_annealer(theta, sampler, mode):
    response = sampler.sample_qubo(theta, num_reads=1)
    l = []
    
    for datum in response.data():
        for key in datum.sample:
            l.append(datum.sample[key])
            pass
        pass

    return l 

def chimera(n):
    G = dnx.chimera_graph(16)
    tmp = nx.to_dict_of_lists(G)
    rows = []
    cols = []
    for i in range(n):
        rows.append(i)
        cols.append(i)
        for j in tmp[i]:
            if(j < n):
                rows.append(i)
                cols.append(j)

    return list(zip(rows, cols))

def send_msg(l, mode):
    if mode == -1:
        for j in range(len(l)):
            if(l[j] == 1):
                msg = ("+" + str(l[j])).encode()
            else:
                l[j] = 1
                msg = ("-" + str(l[j])).encode()
            os.write(1, msg)  # write solution on pipe
            pass
    else:
        for j in range(len(l)):
            msg = (" " + str(l[j])).encode()
            os.write(1, msg)


def main():
    signal.signal(signal.SIGINT, handler)
    mode = sys.argv[1]
    n = sys.stdin.read(6)
    n = int(n.split('\x00', 1)[0])
    
    simulation = sys.stdin.read(2)
    simulation = int(simulation.split('\x00', 1)[0])

    #row = 16 * 8
    #if(n >= row): n_cols = 16
    #else: n_cols = n / 8
    #n_rows = int(n / row)
    #if n % row != 0:
    #    n_rows += 1

    #A = chimera(n_rows, n_cols, n)
    A = chimera(n)

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
        sampler = DWaveSampler()
        sampler = EmbeddingComposite(sampler)
        theta = {}
        i = 0

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
                    theta[(r, c)] = float(x)
                    pass

                i = (i + 1) % 3

            elif(x[0] == "#"):  # compute problem
                i = 0
                l = run_annealer(theta, sampler, mode)
                #l = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

                send_msg(l, mode)

                theta = {}  # clear dictionary
                pass
            
            else:
                end = True
                sys.stderr.write("Closing Python\n")
                pass
            pass


if __name__ == '__main__':
    main()
