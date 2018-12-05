import random
from generator import DistributionGenerator
import config

def generate_populations(num_pops):
    """
    Generate num_pops populations of the type dist_type.
    """
    distGenerator = DistributionGenerator()
    populations = []
    for i in range(num_pops):
        dist_type = random.sample(config.dist_types, 1)[0]
        populations.append(distGenerator.generateDistributions(dist_type, config.MaxDistributionSize))
    return populations

def generate_alphas(num_pops):
    list_of_lists = []
    for i in range(num_pops):
        alphas_list = []
        for j in range(num_pops):
            alphas_list.append(random.random())
        list_of_lists.append(alphas_list)
    return list_of_lists


if __name__ == "__main__":
    print(generate_alphas(5))