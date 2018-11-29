from generator import DistributionGenerator
from sampler import Sampler
from config import *

def multiGT():
    '''
     good toulmin estimator for multi population distriubtions
     returns estimate
    '''

def generateFingerprint(freqList, lists, samples):
    indices = generateIndicies(lists)
    fingerprint = {}
    for index in indices:
        # ilen = len(index)
        if 0 not in index:
            listSet = []
            for i in range(len(index)):
                if index[i] not in freqList[i].keys():
                    fingerprint[index] = 0
                    listSet = []
                    break
                else:
                    listSet.append(set(freqList[i][index[i]]))
            if len(listSet) > 0:
                l = set.intersection(*listSet)
                fingerprint[index] = len(l)
        else:
            zeroList = []
            nonZeroList = []
            for i in range(len(index)):
                if index[i] == 0:
                    zeroList.append(set(uniqueSamples[i]))
                else:
                    if index[i] not in freqList[i]:
                        nonZeroList.append(set())
                    else:
                        nonZeroList.append(set(freqList[i][index[i]]))            
            # import pdb;pdb.set_trace()
            # print 'index: ', index, 'fingerprint', fingerprint
            zeroSet = set.intersection(*zeroList) if len(zeroList) != 0 else set()
            if len(nonZeroList) != 0:
                nonZeroSet = set.intersection(*nonZeroList)
            else:
                nonZeroSet = set()
            fingerprint[index] = len(nonZeroSet.difference(set.intersection(zeroSet, nonZeroSet)))
    return fingerprint

def generateIndicies(lists):
    if any([not l for l in lists]):  
        return
    n = len(lists)
    indexes = [0] * n
    while True:
        yield tuple(lists[i][indexes[i]] for i in xrange(n))  
        for i in xrange(n-1, -1, -1):  
            if indexes[i] < len(lists[i]) - 1:      
                indexes[i] += 1                     
                indexes[i+1:n] = [0] * (n - i - 1)  
                break
        else:
            break

def generateDistributionMappings(distribution):
    countFreqs = {}
    for species in distribution:
        if species not in countFreqs:
            countFreqs[species] = 1
        else:
            countFreqs[species] += 1
    frequencyCounts = {}
    for eachSpecies in countFreqs.keys():
        if countFreqs[eachSpecies] not in frequencyCounts:
            frequencyCounts[countFreqs[eachSpecies]] = [eachSpecies,]
        else:
            l = frequencyCounts[countFreqs[eachSpecies]]
            l.append(eachSpecies)
            frequencyCounts[countFreqs[eachSpecies]] = l
    return max(frequencyCounts, key=int), frequencyCounts

def generateUniqueSamples(distribution):
    unique_set = []
    for val in distribution:
        if val not in unique_set:
            unique_set.append(val)
    return unique_set

if __name__ == "__main__":
    # call the distribution generator to generate populations for various distributions
    # distGenerator = DistributionGenerator()
    # sigmas = distGenerator.generateDistributions()

    # call the sampler for sampling from the above generated distributions
    # sampler = Sampler()
    # samplesX = sampler.generateXSamples()
    # samplesX = sigmas
    # samplesX = [[1,1,2,3,3,3,3,4], [1,2,3,4], [1,3,4,4,4,5]]
    samplesX = [[1,2,3,5,6], [1,2,4,5,5,6,6]]
    freqList = []
    maxKeys = []
    uniqueSamples = []
    for eachDist in samplesX:
        uniqueSamples.append(generateUniqueSamples(eachDist))
        maxKey, distFreq = generateDistributionMappings(eachDist)
        freqList.append(distFreq)
        maxKeys.append([i for i in range(0, maxKey+1)])
    # print maxKeys
    # for the samples generated above, create the multi population fingerprint
    fingerprint = generateFingerprint(freqList, maxKeys, uniqueSamples)
    for i in fingerprint:
        print i, fingerprint[i]
    # implement Good toulmin on the above generated X samples
    # UCap = MultiGT(fingerprint)
    # generate Y samples used to compare with the output of good toulmin
    # samplesY = sampler.generateYsamples()

    #compare and provide some sort of statistic to compare sampleY and MultiGT output

    #generate the joint frequency distribution histogram based on the algorithm given in the paper
    # did not write any skeleton for this part of the code

    #generate some graphs

    print "main file"
