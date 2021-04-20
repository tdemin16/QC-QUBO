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
    return response.first.sample.values()


def run_hybrid(theta):
    sampler = LeapHybridSampler()
    response = sampler.sample_qubo(theta)
    response = response.first.sample.values()

    return response


def decode_solution(tsp_matrix, response, validate=True):
    n = len(tsp_matrix)
    solution = []
    sw = []
    x = []
    last_j = -1
    all = set()
    ins = set()

    for j in response:
        x.append(j)

    for i in range(n):
        solution.append(-1)
        sw.append(0)

    for i in range(n):
        for j in range(n):
            if x[n * i + j] == 1:
                solution[i] = j
                last_j = j
                pass
            pass

        all.add(i)
        if last_j != -1:
            ins.add(last_j)

        last_j = -1

    if validate:
        res = all - ins

        for i in range(n):
            if solution[i] == -1:
                solution[i] = res.pop()
                ins.add(solution[i])

        if len(ins) != len(all):
            for i in range(n):
                if sw[solution[i]] == 0:
                    sw[solution[i]] += 1
                else:
                    solution[i] = res.pop()
                    sw[solution[i]] += 1
        pass

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


def main(n):
    k = 1000

    qubo = dict()
    nodes_array = create_nodes_array(n)
    print(nodes_array, "\n")
    tsp_matrix = get_tsp_matrix(nodes_array)

    constraint_constant = tsp_matrix.max()*len(tsp_matrix)
    cost_constant = 1

    add_cost_objective(tsp_matrix, cost_constant, qubo)
    add_time_constraints(tsp_matrix, constraint_constant, qubo)
    add_position_constraints(tsp_matrix, constraint_constant, qubo)

    start = time()
    solution = solve_tsp(qubo, k)
    end = time()

    solution = decode_solution(tsp_matrix, solution, True)
    cost = calculate_cost(tsp_matrix, solution)

    print("D-Wave solution")
    print("D-Wave:", solution, cost)
    print("Calculation time:", end - start)


def main_hybrid():
    nodes_array = [
        [7.21362447, 4.72776071],
        [7.48093569, 8.15635289],
        [9.72339245, 7.86815924],
        [2.68741026, 2.43794644],
        [6.86692814, 2.86802978],
        [8.11962058, 4.72617735],
        [9.66730978, 1.6601765],
        [0.48265261, 9.97703727],
        [4.85291697, 3.57640542],
        [4.64010762, 0.77011666]
    ]

    qubo = dict()
    tsp_matrix = get_tsp_matrix(nodes_array)
    constraint_constant = tsp_matrix.max()*len(tsp_matrix)
    cost_constant = 1

    add_cost_objective(tsp_matrix, cost_constant, qubo)
    add_time_constraints(tsp_matrix, constraint_constant, qubo)
    add_position_constraints(tsp_matrix, constraint_constant, qubo)

    start = time()
    solution = run_hybrid(qubo)
    end = time()

    solution = decode_solution(tsp_matrix, solution, True)
    cost = calculate_cost(tsp_matrix, solution)

    print("Hybrid solution")
    print("Hybrid: ", solution, cost)
    print("Calculation time:", end - start)



if __name__ == "__main__":
    main_hybrid()
