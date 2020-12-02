import time
from dwave.system.samplers import DWaveSampler
from neal import SimulatedAnnealingSampler
from dwave.system.composites import EmbeddingComposite

def dec(n):
    num = n
    dec_value = 0

    # Initializing base
    # value to 1, which is 2 ^ 0
    base1 = 1

    len1 = len(num)
    for i in range(len1 -1, -1, -1):
        if (num[i] == '1'):
            dec_value += base1
            base1 = base1 * 2
            pass
        pass

    return dec_value

useQpu = False

if(useQpu):
    sampler = DWaveSampler()
    sampler = EmbeddingComposite(sampler)
else:
    sampler = SimulatedAnnealingSampler()

bqm = {}
distrib = {} #don't think I need

bits = 4

for i in range(0, bits):
    bqm[(i,i)] = 0
    pass

start = time.time()
response = sampler.sample_qubo(bqm, num_reads=1)
total = time.time() - start

num = ""
for datum in response.data():
    for key in datum.sample:
        num += str(datum.sample[key])

print(dec(num))
