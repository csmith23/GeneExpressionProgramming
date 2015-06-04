__author__ = 'Coleman'

import random
import math

class Environment:
    def __init__(self, populationSize, mutationRate=0, inversionRate=0, ISTranspositionRate=0, RISTranspositionRate=0, geneTranspositionRate=0, onePointRecombinationRate=0, twoPointRecombinationRate=0, geneRecombinationRate=0):
        self.population = []
        self.populationSize = populationSize
        self.mutationRate = mutationRate
        self.inversionRate = inversionRate
        self.ISTranspositionRate = ISTranspositionRate
        self.RISTranspositionRate = RISTranspositionRate
        self.geneTranspositionRate = geneTranspositionRate
        self.onePointRecombinationRate = onePointRecombinationRate
        self.twoPointRecombinationRate = twoPointRecombinationRate
        self.geneRecombinationRate = geneRecombinationRate

    def init(self, headLengths, homeoticLengths, genomes):
        for i in range(self.populationSize):
            self.population.append(Chromosome(headLengths, homeoticLengths, genomes))

    def run(self, fitness):
        pass
        # get fitness
        # select most fit
        # mutate new generation

    def select(self):
        pass
        # roulette wheel selection, with most fit guaranteed

    def reproduce(self):
        pass

    def alter(self):
        pass

class Chromosome:
    def __init__(self, headLengths, homeoticLengths, genomes):
        self.fitness = 0
        self.symbols = [genomes[i].symbols(len(headLengths), i) for i in range(len(genomes))]
        self.genomes = genomes

        terminals = [genomes[i].terminals for i in range(len(headLengths))]
        homeoticTerminals = [list(range(len(headLengths))) for i in range(len(homeoticLengths))]
        self.terminals = terminals + homeoticTerminals

        heads = [genomes[i].head(headLengths[i]) for i in range(len(headLengths))]
        tails = [genomes[i].tail(headLengths[i]) for i in range(len(headLengths))]
        homeotics = [genomes[i + len(headLengths)].homeotic(homeoticLengths[i], len(headLengths)) for i in range(len(homeoticLengths))]
        homeoticTails = [genomes[i + len(headLengths)].homeoticTail(homeoticLengths[i], len(headLengths)) for i in range(len(homeoticLengths))]
        self.headRanges = []
        self.tailRanges = []
        self.homeoticRanges = []
        self.sequence = []
        for i in range(len(heads)):
            self.headRanges.append((len(self.sequence), len(self.sequence) + len(heads[i])))
            self.tailRanges.append((len(self.sequence) + len(heads[i]), len(self.sequence) + len(heads[i]) + len(tails[i])))
            self.homeoticRanges.append(False)
            self.sequence = self.sequence + heads[i]
            self.sequence = self.sequence + tails[i]

        for i in range(len(homeotics)):
            self.headRanges.append((len(self.sequence), len(self.sequence) + len(homeotics[i])))
            self.tailRanges.append((len(self.sequence) + len(homeotics[i]), len(self.sequence) + len(homeotics[i]) + len(homeoticTails[i])))
            self.homeoticRanges.append(True)
            self.sequence = self.sequence + homeotics[i]
            self.sequence = self.sequence + homeoticTails[i]

        self.head = [False for i in range(len(self.sequence))]
        self.tail = [True for i in range(len(self.sequence))]
        self.homeotic = [False for i in range(len(self.sequence))]
        for headRange in self.headRanges:
            for i in range(headRange[0], headRange[1]):
                self.head[i] = True
                self.tail[i] = False

        for index in range(len(self.homeoticRanges)):
            if self.homeoticRanges[index]:
                for i in range(self.headRanges[index][0], self.tailRanges[index][1]):
                    self.homeotic[i] = True

    def eval(self):
        pass

    def mutation(self, times):
        indices = list(range(len(self.sequence)))
        mutationIndices = [indices.pop(random.randint(0, len(indices) - 1)) for i in range(times)]
        print(end="\n\n\n")
        print(mutationIndices)
        for index in mutationIndices:
            if self.head[index]:
                self.sequence[index] = False

    def inversion(self):
        pass

    def ISTransposition(self):
        pass

    def RISTransposition(self):
        pass

    def geneTransposition(self):
        pass

    def onePointRecombination(self):
        pass

    def twoPointRecombination(self):
        pass

    def geneRecombination(self):
        pass

class Genome:
    def __init__(self):
        self.functions = {}
        self.arities = {}
        self.terminals = []
        self.arity = 0

    def addFunction(self, symbol, function, arity):
        self.functions[symbol] = function
        self.arities[symbol] = arity
        if arity > self.arity:
            self.arity = arity

    def addTerminal(self, symbol):
        self.terminals.append(symbol)

    def head(self, length):
        head = [random.choice(list(self.functions.keys()) + self.terminals) for i in range(length)]
        head[0] = random.choice(list(self.functions.keys()))
        return head

    def tail(self, length):
        tail = [random.choice(self.terminals) for i in range(length * (self.arity - 1) + 1)]
        return tail

    def homeotic(self, length, numGenes):
        homeotic = [random.choice(list(self.functions.keys()) + list(range(numGenes))) for i in range(length)]
        homeotic[0] = random.choice(list(self.functions.keys()))
        return homeotic

    def homeoticTail(self, length, numGenes):
        tail = [random.choice(list(range(numGenes))) for i in range(length * (self.arity - 1) + 1)]
        return tail

    def symbols(self, numGenes, index):
        if index > numGenes - 1:
            return list(self.functions.keys()) + list(range(numGenes))

        return list(self.functions.keys()) + self.terminals


    def _add(inputs):
        if inputs[0] == "Error" or inputs[1] == "Error":
            return "Error"
        return inputs[0] + inputs[1]
    def _sub(inputs):
        if inputs[0] == "Error" or inputs[1] == "Error":
            return "Error"
        return inputs[0] - inputs[1]
    def _mul(inputs):
        if inputs[-1] == "Error" or inputs[-2] == "Error":
            return "Error"
        return inputs[0] * inputs[1]
    def _div(inputs):
        if inputs[0] == "Error" or b == "Error":
            return 1
        elif inputs[1] != 0:
            return inputs[0] / inputs[1]
        else:
            return "Error"

    def _NOT(inputs):
        return not inputs[0]
    def _AND(inputs):
        return inputs[0] and inputs[1]
    def _OR(inputs):
        return inputs[0] or inputs[1]
    def _XOR(inputs):
        return inputs[0] is not inputs[1]
    def _XNOR(inputs):
        return inputs[0] is inputs[1]

    def _sqr(inputs):
        return inputs[0] ** 2
    def _sqrt(inputs):
        return inputs[0] ** 0.5

    def _exp(inputs):
        return 2 ** inputs[0]

    def _log(inputs):
        return math.log(inputs[0], 2)

    ARITHMETIC_SET = ({"+": _add, "-": _sub, "*": _mul, "/": _div}, {"+": 2, "-": 2, "*": 2, "/": 2})
    BOOLEAN_SET = ({"N": _NOT, "A": _AND, "O": _OR}, {"N": 2, "A": 2, "O": 2})
    REEDMULLEN_SET = ({"A": _AND, "X": _XOR, "N": _NOT}, {"A": 2, "X": 2, "N": 2})
    QUADRATIC_SET = ({"E": _sqr, "L": _sqrt}, {"E": 1, "L": 1})
    EXPONENTIAL_SET = ({"E": _exp, "L": _log}, {"E": 1, "L": 1})

class Node:
    def __init__(self, symbol):
        self.children = []
        self.symbol = symbol

    def eval(self, genome, terminalInputs, homeoticInputs):
        inputs = []
        for child in self.children:
            inputs.append(child.eval(genome, terminalInputs, homeoticInputs))

        if self.symbol in list(genome.functions.keys()):
            return genome.functions[self.symbol](inputs)
        elif isinstance(self.symbol, int):
            return homeoticInputs[self.symbol]

        return terminalInputs[self.symbol]


genome = Genome()
genome.addFunction("A", None, 2)
genome.addFunction("O", None, 2)
genome.addFunction("N", None, 1)
genome.addTerminal("a")
genome.addTerminal("b")
chromosome = Chromosome([3, 3, 3], [3], [genome for i in range(4)])
print(chromosome.headRanges)
print(chromosome.sequence)
print(chromosome.homeoticRanges)
print(chromosome.tailRanges)
print(chromosome.symbols)
print(chromosome.head)
print(chromosome.tail)
print(chromosome.homeotic)