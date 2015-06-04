__author__ = 'Coleman'

import abc
import random

class Environment(metaclasss = abc.Meta):
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

class Chromosome():
    def __init__(self):
        self.genes = []
        self.link = ""

    def initRand(self, numGenes, headLength, genome, link):
        self.genes = [Gene(genome).initRand(headLength) for i in range(numGenes)]

    def replicate(self):
        newGenes = []
        for gene in self.genes:
            newGenes.append(gene.replicate())

        return newGenes

class Gene():
    def __init__(self, genome):
        self.genome = genome
        self.head = []
        self.tail = []

    def initRand(self, length):
        arity = self.genome.functions["arity"]
        functions = self.genome.functions.keys()
        functions.remove("arity")
        terminals = self.genome.terminals
        self.head = [random.choice(functions + terminals) for i in range(length)]
        self.tail = [random.choice(terminals) for i in range(length * (arity - 1) + 1)]
        return self

    def eval(self):
        '''only for use if the eval algorithm works, and is implemented, otherwise use the inefficient and error-prone node system'''
        pass

    def replicate(self):
        newGene = Gene(self.genome)
        newGene.head = self.head[:]
        self.tail = self.tail[:]
        return newGene

class Genome():
    def _add(inputs):
        return sum(inputs)
    def _sub(inputs):
        return inputs[0] - inputs[1]
    def _mul(inputs):
        return inputs[0] * inputs[1]
    def _div(inputs):
        return inputs[0] / inputs[1]

    def _NOT(inputs):
        return not inputs[0]
    def _AND(inputs):
        return inputs[0] and inputs[1]
    def _OR(inputs):
        return inputs[0] or inputs[1]
    def _XOR(inputs):
        return inputs[0] != inputs[1]
    def _XNOR(inputs):
        return inputs[0] == inputs[1]

    def _exp(inputs):
        return inputs[0] ** 2
    def _sqrt(inputs):
        return inputs[0] ** 0.5

    ARITHMETIC_SET = {"+": (_add, 2), "-": (_sub, 2), "*": (_mul, 2), "/": (_div, 2), "arity": 2}
    BOOLEAN_SET = {"N": (_NOT, 1), "A": (_AND, 2), "O": (_OR, 2), "arity": 2}
    REEDMULLEN_SET = {"A": (_AND, 2), "X": (_XOR, 2), "N": (_NOT, 1), "arity": 2}
    EXPONENTIAL_SET = {"E": (_exp, 1), "Q": (_sqrt, 1), "arity": 1}

    def __init__(self):
        self.functions = {"arity": 0}
        self.terminals = []