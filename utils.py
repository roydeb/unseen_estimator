import random
from generator import DistributionGenerator
import config
import math as mt
import pandas as pd

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

def generate_prob_distributions(pops):
    pop_dists = []
    for pop in pops:
        s = pd.Series(pop)
        pop_dists.append(list((s.groupby(s).transform('count') / len(s)).values))
    return pop_dists


def binomial(p, n, i):
    log = mt.log
    exp = mt.exp
    output = n*log(1-p)
    for i in range(0,i):
        output += log(n-i) + log(p/(1-p)) - log(i+1)
    return exp(output)


# if __name__ == "__main__":
#     print(generate_alphas(5))