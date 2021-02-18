import dimod
import hybrid
import random
from dwave.system.samplers import LeapHybridSampler
import numpy as np
import time
from math import sqrt


def run_annealer_hybrid(theta):
    sampler = LeapHybridSampler()
    response = sampler.sample_qubo(theta)
    response = response.first.sample.values()

    return np.atleast_2d(list(response)).T


def to_NPP(nums):
    theta = [[0 for col in range(len(nums))] for row in range(len(nums))]
    c = sum(nums)

    for i in range(len(nums)):
        for j in range(len(nums)):
            if i != j:
                theta[i][j] = nums[i] * nums[j]
            else:
                theta[i][i] = nums[i] * (nums[i] - c)

    return np.array(theta)


def get_QAP(namefile, lam):
    fo = open(namefile, "r")
    n = int(fo.read(1))

    i = 0
    count = 0

    f = [[0 for col in range(n)] for row in range(n)]
    d = [[0 for col in range(n)] for row in range(n)]
    line = fo.read()
    for word in line.split():
        if count < 9:
            f[int(i/3)][int(i % 3)] = int(word)
        else:
            d[int(i/3)][int(i % 3)] = int(word)

        i = i + 1
        count += 1
        if(count == 9):
            i = 0
        pass

    fo.close()

    Q = np.kron(f, d)
    m = Q.max()
    pen = lam * m
    Q = np.einsum("ij,kl->ikjl", f, d).astype(np.float)
    i = range(len(Q))

    Q[i, :, i, :] += pen
    Q[:, i, :, i] += pen
    Q[i, i, i, i] -= 4 * pen

    Q = Q.reshape(n**2, n**2)

    return Q, m, pen


def matrix_to_dict(theta):
    n = len(theta)
    d = dict()
    for i in range(n):
        for j in range(n):
            d[i, j] = theta[i][j]

    return d


def gen(n, ran):
    l = []
    for i in range(int(n)):
        l.append(random.randint(0, ran))

    return l


def fQ(theta, sol):
    return ((np.atleast_2d(sol).T).dot(theta)).dot(sol)


namefile = "./test/test_1.txt"
theta, m, pen = get_QAP(namefile, 1.25)

start = time.time()
print("Start hybrid")
sol = run_annealer_hybrid(matrix_to_dict(theta))
end = time.time()

print(str(end - start) + "s")
print("max coeff=" + str(m))
print("y=" + str(pen * sqrt(len(theta)) * 2 + fQ(theta, sol)))
print("pen=" + str(pen))
print("fQ=" + str(fQ(theta, sol)))