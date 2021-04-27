#!/usr/bin/python3
#
#
#
#       THE CODE IN THIS PYTHON SCRIPT IS NOT FROM US
#           WE ARRANGED AN ALREADY EXISTING SCRIPT
#   FROM: https://github.com/BOHRTECHNOLOGY/quantum_tsp
#
#
#


import numpy as np
import sys
# Library to interact with the QPU
from dwave.system.samplers import DWaveSampler
# Library to embed our problem onto the QPU physical graph
from dwave.system.composites import EmbeddingComposite
from dwave.system.samplers import LeapHybridSampler
import random
from time import time


def add_cost_objective(distance_matrix, cost_constant, qubo_dict):
    n = len(distance_matrix)
    for t in range(n):
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                qubit_a = t * n + i
                qubit_b = (t + 1) % n * n + j
                qubo_dict[(qubit_a, qubit_b)] = cost_constant * \
                    distance_matrix[i][j]


def add_time_constraints(distance_matrix, constraint_constant, qubo_dict):
    n = len(distance_matrix)
    for t in range(n):
        for i in range(n):
            qubit_a = t * n + i
            if (qubit_a, qubit_a) not in qubo_dict.keys():
                qubo_dict[(qubit_a, qubit_a)] = -constraint_constant
            else:
                qubo_dict[(qubit_a, qubit_a)] += -constraint_constant
            for j in range(n):
                qubit_b = t * n + j
                if i != j:
                    qubo_dict[(qubit_a, qubit_b)] = 2 * constraint_constant


def add_position_constraints(distance_matrix, constraint_constant, qubo_dict):
    n = len(distance_matrix)
    for i in range(n):
        for t1 in range(n):
            qubit_a = t1 * n + i
            if (qubit_a, qubit_a) not in qubo_dict.keys():
                qubo_dict[(qubit_a, qubit_a)] = -constraint_constant
            else:
                qubo_dict[(qubit_a, qubit_a)] += -constraint_constant
            for t2 in range(n):
                qubit_b = t2 * n + i
                if t1 != t2:
                    qubo_dict[(qubit_a, qubit_b)] = 2 * constraint_constant


def solve_tsp(qubo_dict, k):
    response = EmbeddingComposite(DWaveSampler()).sample_qubo(
        qubo_dict, chain_strength=800, num_reads=k)
    return list(response.first.sample.values())


def run_hybrid(theta):
    sampler = LeapHybridSampler()
    response = sampler.sample_qubo(theta)
    response = response.first.sample.values()

    return list(response)


def advance(iter, rnd):
    iterator = next(iter)
    while random.random() > rnd:
        iterator = next(iter, iterator)
    return iterator


def decode_solution(response, validate):
    
    n = int(np.sqrt(len(response)))
    solution = np.array(n)
    raw = dict()
    for i in range(n):
        raw[i] = list()
    keep = list()
    all_ = list()
    diff = list()
    indexes = list()

    if not validate:
        solution = list()
        for i in range(n):
            for j in range(n):
                if(response[n*i + j] == 1):
                    last = j
            if (last != -1):
                solution.append(last)
            last = -1
    else:
        solution = np.array([-1 for i in range(n)])
        for i in range(n):
            for j in range(n):
                if (response[n*i +j] == 1):
                    raw[i].append(j)
        
        for i in range(n):
            if len(raw[i]) == 1:
                keep.append(raw[i][0])
                solution[i] = raw[i][0]
            all_.append(i)            

        for i in range(n):
            if len(raw[i]) > 1:
                for it in raw[i]:
                    if it not in keep: 
                        diff.append(it)
                
                if len(diff) > 0:
                    it = advance(iter(diff), random.random() % len(diff))
                    solution[i] = it
                    keep.append(it)
                    diff.clear()

        for i in range(n):
            for j in range(n):
                if solution[j] == i:
                    indexes.append(j) 

            if len(indexes) > 1:
                random.shuffle(indexes)
                index = indexes[0]
                for it in indexes:
                    if it == index: 
                        solution[it] = i 
                    else:
                        solution[it] = -1

                keep.append(i)

            indexes.clear()

        for it in all_:
            if it not in keep: 
                diff.append(it)
            
        for i in range(n):
            if solution[i] == -1 and len(diff) != 0:
                it = advance(iter(diff), random.random() % len(diff))
                solution[i] = it
                diff.remove(it)

    return solution


# def decode_solution(distance_matrix, response):
#     best_solution = np.inf
#     n = len(distance_matrix)
#     min_energy = response.record[0].energy
#     for record in response.record:
#         sample = record[0]
#         solution_binary = [node for node in sample]
#         solution = binary_state_to_points_order(solution_binary)
#         if record.energy <= min_energy:
#             best_solution = solution

#     return best_solution


def calculate_cost(distance_matrix, solution):
    cost = 0
    for i in range(len(solution)):
        a = i % len(solution)
        b = (i + 1) % len(solution)
        cost += distance_matrix[solution[a]][solution[b]]

    return cost


def binary_state_to_points_order(binary_state):
    points_order = []
    number_of_points = int(np.sqrt(len(binary_state)))
    for p in range(number_of_points):
        for j in range(number_of_points):
            if binary_state[(number_of_points) * p + j] == 1:
                points_order.append(j)
    return points_order


def create_nodes_array(N):
    nodes_list = []
    for i in range(N):
        nodes_list.append(np.random.rand(2) * 10)
    return np.array(nodes_list)


def get_tsp_matrix(nodes_array):
    n = len(nodes_array)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            matrix[i][j] = distance(nodes_array[i], nodes_array[j])
            matrix[j][i] = matrix[i][j]
    return matrix


def distance(point_A, point_B):
    return np.sqrt((point_A[0] - point_B[0])**2 + (point_A[1] - point_B[1])**2)


def main(n, old):
    k = 1000

    qubo = dict()
    if old == 0:
        nodes_array = [
 [9.94842154, 7.98192104],
 [8.60849611, 8.00248268],
 [7.7510664 , 8.92952639],
 [2.48323566, 7.71263975],
 [0.79106925, 4.7339348 ],
 [5.43068116, 8.76787367],
 [2.01200685, 7.57374064],
 [8.24493747, 4.47630442],
 [6.12997546, 7.6483757 ],
 [5.42759893, 6.52504111],
 [1.19438987, 6.82360516],
 [5.30828759, 3.38905631]]
    else:
        nodes_array = create_nodes_array(n)


    print(nodes_array, "\n")
    tsp_matrix = get_tsp_matrix(nodes_array)

    constraint_constant = tsp_matrix.max()*len(tsp_matrix)
    cost_constant = 1

    add_cost_objective(tsp_matrix, cost_constant, qubo)
    add_time_constraints(tsp_matrix, constraint_constant, qubo)
    add_position_constraints(tsp_matrix, constraint_constant, qubo)

    start = time()
    solution = run_hybrid(qubo)
    end = time()

    solution = decode_solution(solution, True)
    cost = calculate_cost(tsp_matrix, solution)

    print("Hybrid solution")
    print("Hybrid: ", solution, cost)
    print("Calculation time:", end - start)

    
    start = time()
    solution = solve_tsp(qubo, k)
    end = time()

    solution = decode_solution(solution, True)
    cost = calculate_cost(tsp_matrix, solution)

    print("D-Wave solution")
    print("D-Wave:", solution, cost)
    print("Calculation time:", end - start)


if __name__ == "__main__":
    n = int(input("Insert n: "))
    old = int(input("New dataset?(0/1):"))
    main(n, old)
