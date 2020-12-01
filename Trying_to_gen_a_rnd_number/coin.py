import time
from dwave.system.samplers import DWaveSampler
from neal import SimulatedAnnealingSampler
from dwave.system.composites import EmbeddingComposite

useQpu = False

if(useQpu):
    sampler = DWaveSampler()
    sampler = EmbeddingComposite(sampler)
else:
    sampler = SimulatedAnnealingSampler()

bqm = {}
distrib = {} #don't think I need
