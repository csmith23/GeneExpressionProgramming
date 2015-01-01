__author__ = 'Coleman'

import abc
import random
import copy

from GEP.Genome import Genome
from GEP.Chromosome import *

class Environment(metaclass = abc.ABCMeta):
    def __init__(self):
        self.chromosomes = [] # pair of chromosome and the its fitness value
        self.inputsOutputs = [] # pair of list of inputs and expected output value
        '''[([0, 0, 0], 0),
            ([0, 0, 1], 0),
            ([0, 1, 0], 0),
            ([0, 1, 1], 1),
            ([1, 0, 0], 0),
            ([1, 0, 1], 1),
            ([1, 1, 0], 1),
            ([1, 1, 1], 1),]'''

        self.generations = 50
        self.populationSize = 20
        self.headLength = 7
        self.genome = Genome()
        self.numGenes = 3
        self.linkingGene = Gene()

        self.mutationRate = 0.044
        self.inversionRate = 0.1
        self.ISTranspositionRate = 0.1
        self.RISTranspositionRate = 0.1
        self.GeneTranspositionRate = 0.1
        self.OnePointRecombination = 0.4
        self.TwoPointRecombinationRate = 0.2
        self.GeneRecombinationRate = 0.1

    def setup(self):
        '''all instance attributes in the second group must be initialize beforehand'''

        self.chromosomes = [Chromosome(self.genome).initRand(self.numGenes, self.headLength) for i in range(self.populationSize)]
        for chromosome in self.chromosomes:
            chromosome.link = self.linkingGene

    def run(self):
        '''targetOutputs = [self.inputsOutputs[i][1] for i in range(len(self.inputsOutputs))]
           inputs = [self.inputsOutputs[i][0] for i in range(len(self.inputsOutputs))]
           if [a, b, c] is the list of terminals in self.genome, inputs [1, 2, 3] corresponds to a = 1, b = 2, c = 3'''
        targetOutputs = [self.inputsOutputs[i][1] for i in range(len(self.inputsOutputs))]
        outputs = []
        inputs = [self.inputsOutputs[i][0] for i in range(len(self.inputsOutputs))]
        keys = self.genome.terminals

        for i in range(len(self.chromosomes)):
            chromosome = self.chromosome[i]

            chromosome[0].Tree()

            for input in inputs:
                if len(input) != len(keys):
                    return

                terminalDict = dict(zip(keys, input))
                outputs.append = chromosome[0].tree.step(terminalDict)

            self.chromosomes[i] = (chromosome[0], self.fitness(outputs, targetOutputs))

        self.select()
        self.modify()

    @abc.abstractmethod
    def fitness(self, actualOutputs, targetOutputs):
        '''returns a fitness value based on error, called on each element in self.chromosomes by run'''
        pass

    def select(self):
        '''roulette wheel based selection for population size - number of chromosomes with best fitness of generation
           the chromosomes with best fitness of generation are automatically selected to be continued'''
        newChromosomes = []
        wheel = []
        for i in range(self.chromosomes):
            for j in range(self.chromosomes[i][1]):
                wheel.append(i)
        random.shuffle(wheel)

        eliteFitness = 0
        for chromosome in self.chromosomes:
            if chromosome[1] > eliteFitness:
                eliteFitness = chromosome[1]

        numElite = 0
        for i in range(self.chromosomes):
            if self.chromosome[i][1] == eliteFitness:
                newChromosomes.append(i)
                numElite += 1

        for i in range(self.populationSize - numElite):
            newChromosomes.append(random.choice(wheel))

        self.reproduce(newChromosomes)

    def reproduce(self, numDaughters):
        '''given a list of chromosomes and the number of times for the to be duplicated,
           sets new generation of randomly order list of all the duplicated chromosomes
           the tuples in self.chromosome are create here, not in select'''
        pass

    def modify(self):
        '''calls the genetic modifiers with their respective probabilities'''
        for chromosome in self.chromosomes:
            chromosome[0].mutation(self.mutationRate)
            chromosome[0].inversion(self.inversionRate)
            chromosome[0].ISTransposition(self.ISTranspositionRate)
            chromosome[0].RISTransposition(self.RISTranspositionRate)
            chromosome[0].geneTransposition(self.GeneTranspositionRate)
            chromosome[0].onePointRecombination(self.OnePointRecombination)
            chromosome[0].twoPointRecombination(self.TwoPointRecombinationRate)
            chromosome[0].geneRecombination(self.GeneRecombinationRate)
