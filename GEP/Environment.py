__author__ = 'Coleman'

import abc
from Chromosome import Chromosome
from Genome import Genome
from Gene import Gene

class Environment(metaclass=abc.meta):
    def __init__(self):
        self.chromosomes = []

        self.populationSize = 0
        self.headLength = 0
        self.link = ""

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
        pass

    def printChromosomes(self):
        pass