import random
from generator import DistributionGenerator
import config

def generate_populations(num_pops, dist_type):
    """
    Generate num_pops populations of the type dist_type.
    """
    distGenerator = DistributionGenerator()
    populations = []
    for i in range(num_pops):
        populations.append(distGenerator.generateDistributions(dist_type, config.MaxDistributionSize))
    return populations