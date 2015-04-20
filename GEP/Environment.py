__author__ = 'Coleman'

import abc
from Chromosome import Chromosome
from Genome import Genome
from Gene import Gene

class Environment(metaclass=abc.meta):
    def __init__(self):
        self.population = []

        self.populationSize = 0
        self.headLength = 0
        self.homeotic = None

        self.homeoticRate = 0 # multiplied by the other rates for every homeotic gene
        self.mutationRate = 0
        self.inversionRate = 0
        self.ISTranspositionRate = 0
        self.RISTranspositionRate = 0
        self.geneTranspositionRate = 0
        self.onePointRecombinationRate = 0
        self.twoPointRecombinationRate = 0
        self.geneTranspositionRate = 0

    def run(self):
        pass

    @abc.abstractmethod
    def fitness(self):
        pass

    def select(self):
        pass

    def replicate(self, indicies):
        pass

    def modify(self):
        for chromosome in self.population:
            pass

    def printChromosomes(self):
        pass