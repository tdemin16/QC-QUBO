import dimod
import os
import sys
import dwave_networkx as dnx
import networkx as nx
import scipy
import signal
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite


# Handler sigint
def handler(signum, frame):
    exit(0)

# Given a theta dictionary and a sampler, computes the global minimum
def run_annealer(theta, sampler, k):
    # Run the annealer 4 times with theta matrix
    response = sampler.sample_qubo(theta, num_reads=k)
    
    # Samples are orderes from lowest energy to highest -> fist sample has lowest energy
    response = response.first.sample.values()
    sys.stderr.write(str(response) + "\n")

    return response

# Given n number of nodes, generates a chimera graph of n nodes
def chimera(n):
    G = dnx.chimera_graph(16)
    tmp = nx.to_dict_of_lists(G)
    rows = []
    cols = []

    for i in range(n):
        
        # Append the diagonal
        rows.append(i)
        cols.append(i)

        # Append each edge within the predefined range
        for j in tmp[i]:
            if (j < n):
                rows.append(i)
                cols.append(j)

    return list(zip(rows, cols))

# Given n number of nodes, generates a pegasus graph of n nodes
def pegasus(n):
    G = dnx.pegasus_graph(16, fabric_only=False)
    tmp = nx.to_numpy_matrix(G)
    
    graph = []
           
    for i in range(n):
        graph.append((i, i))

        for j in range(n):
            if tmp.item(i,j):
                graph.append((i, j))

    return graph

# Given a list l and a bool value mode, send l back to C++
def send_msg(l):
    for j in l:
        msg = (" " + str(j)).encode()
        os.write(1, msg)


def send_topology(active_nodes, active_edges, n):
    len_n = len(str(n)) + 5
    
    # send active nodes
    for node in active_nodes:
        msg = str(node)
        msg = ('0' * (len_n - len(msg)) + msg).encode() 
        os.write(1, msg)
        pass

    # send active edges
    for r, c in active_edges:
        msg = str(r)
        msg = ('0' * (len_n - len(msg)) + msg).encode() 
        os.write(1, msg)

        msg = str(c)
        msg = ('0' * (len_n - len(msg)) + msg).encode()
        os.write(1, msg)
        pass       

    # End of transmission
    msg = ("#" * len_n).encode()
    os.write(1, msg)
    os.write(1, msg)


def main():
    signal.signal(signal.SIGINT, handler) # Add signal handling
    
    k = int(sys.argv[1]) # read from argv
    
    n = sys.stdin.read(10)                # Read problem's dimension from stdin (pipe)
    n = int(n.split('\x00', 1)[0])        # Decode the dimension
    
    simulation = sys.stdin.read(2)                   # Read type of run from stdin (could be a simulation or not)
    simulation = int(simulation.split('\x00', 1)[0]) # Decode the type of run
    

    if (simulation != 0):
        active_nodes = [i for i in range(n)]
        active_edges = pegasus(n) # Generate a pegasus graph given a a dimension n
    else:
        sampler = DWaveSampler()
        active_nodes = sampler.nodelist[:n]
        active_edges = []
        for edge in sampler.edgelist:
            if edge[0] <= active_nodes[-1] and edge[1] <= active_nodes[-1]: 
                active_edges.append((edge[0], edge[1]))
                active_edges.append((edge[1], edge[0]))
        
        for node in active_nodes:
            active_edges.append((node, node))

    # Send pegasus topology to C++
    send_topology(active_nodes, active_edges, n)


    if(simulation == 0): # If this run must use the annealer then

        theta = {}
        i = 0

        end = False
        while end == False:
            x = sys.stdin.read(100)
            if(x[0] != "#" and x[0:3] != "END"):  # retrieving problem from pipe
                x = x.split('\x00', 1)[0]
                if i == 0:
                    r = int(x) # row
                elif i == 1:
                    c = int(x) # column
                else:
                    theta[(r, c)] = float(x) # val
                    pass

                i = (i + 1) % 3

            elif(x[0] == "#"):  # if '#' is received, compute the minimum
                i = 0
                l = run_annealer(theta, sampler, k) # run the annealer with theta and sampler

                send_msg(l) # Send the result

                theta = {}  # clear dictionary
            
            else:
                end = True
                sys.stderr.write("Closing Python\n")
                pass
            pass


if __name__ == '__main__':
    main()
