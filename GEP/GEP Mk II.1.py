__author__ = 'Coleman'

import random
import math
import sys

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
        while True:
            # calculate fitness
            self.reproduce(self.select())
            self.alter()
            self.print()

    def select(self):
        wheel = []
        for index in range(self.populationSize):
            for i in range(chromosome.fitness):
                wheel.append(index)

        indices = []
        for i in range(self.populationSize):
            index = random.randrange(len(wheel))
            indices.append(wheel[index])

        return indices

    def reproduce(self, indices):
        newPopulation = []
        for index in indices:
            newPopulation.append(self.population[index].copy())

        self.population = newPopulation

    def alter(self):
        for chromosome in self.population:
            chromosome.mutation(self.mutationRate)

        random.shuffle(self.population)
        for chromosome in self.population[:int(round(self.populationSize * self.inversionRate))]:
            chromosome.inversion()

        random.shuffle(self.population)
        for chromosome in self.population[:int(round(self.populationSize * self.ISTranspositionRate))]:
            chromosome.ISTransposition()

        random.shuffle(self.population)
        for chromosome in self.population[:int(round(self.populationSize * self.RISTranspositionRate))]:
            chromosome.RISTranspoisition()

        random.shuffle(self.population)
        for chromosome in self.population[:int(round(self.populationSize * self.geneTranspositionRate))]:
            chromosome.geneTransposition()

        random.shuffle(self.population)
        for chromosome in self.population[:int(round(self.populationSize * self.onePointRecombinationRate))]:
            while True:
                other = random.choice(self.population)
                if random is not chromosome:
                    chromosome.onePointRecombination()
                    break

        random.shuffle(self.population)
        for chromosome in self.population[:int(round(self.populationSize * self.twoPointRecombinationRate))]:
            while True:
                other = random.choice(self.population)
                if random is not chromosome:
                    chromosome.twoPointRecombination()
                    break

        random.shuffle(self.population)
        for chromosome in self.population[:int(round(self.populationSize * self.geneRecombinationRate))]:
            while True:
                other = random.choice(self.population)
                if random is not chromosome:
                    chromosome.genePointRecombination()
                    break

    def print(self):
        pass

class Chromosome:
    def __init__(self, headLengths, homeoticLengths, genomes, gen=True):
        self.fitness = 0
        # self.terminals = None
        self.symbols = None
        self.genomes = None

        self.headRanges = None
        self.tailRanges = None
        self.homeoticRanges = None
        self.sequence = None

        self.head = None
        if gen:
            # terminals = [genomes[0].terminals for i in range(len(headLengths))]
            # homeoticTerminals = [list(range(len(headLengths))) for i in range(len(homeoticLengths))]
            # self.terminals = terminals + homeoticTerminals
            self.symbols = [genomes[0].symbols(len(headLengths), False)]
            self.symbols.append(genomes[1].symbols(len(headLengths), True))
            self.genomes = genomes

            heads = [genomes[0].head(headLengths[i]) for i in range(len(headLengths))]
            tails = [genomes[0].tail(headLengths[i]) for i in range(len(headLengths))]
            homeotics = [genomes[1].homeotic(homeoticLengths[i], len(headLengths)) for i in range(len(homeoticLengths))]
            homeoticTails = [genomes[1].homeoticTail(homeoticLengths[i], len(headLengths)) for i in range(len(homeoticLengths))]
            self.headRanges = []
            self.tailRanges = []
            self.homeoticRanges = []
            self.sequence = []
            for i in range(len(heads)):
                self.headRanges.append((len(self.sequence), len(self.sequence) + len(heads[i])))
                self.tailRanges.append((len(self.sequence) + len(heads[i]), len(self.sequence) + len(heads[i]) + len(tails[i])))
                self.homeoticRanges.append(False)
                self.sequence += heads[i]
                self.sequence += tails[i]

            for i in range(len(homeotics)):
                self.headRanges.append((len(self.sequence), len(self.sequence) + len(homeotics[i])))
                self.tailRanges.append((len(self.sequence) + len(homeotics[i]), len(self.sequence) + len(homeotics[i]) + len(homeoticTails[i])))
                self.homeoticRanges.append(True)
                self.sequence += homeotics[i]
                self.sequence += homeoticTails[i]

            self.head = [False for i in range(len(self.sequence))]
            for headRange in self.headRanges:
                for i in range(headRange[0], headRange[1]):
                    self.head[i] = True

    def eval(self, terminalInputs):
        genes = []
        for headRange, tailRange, homeoticRange in zip(self.headRanges, self.tailRanges, self.homeoticRanges):
            genes.append((self.sequence[headRange[0]:tailRange[1]], homeoticRange))

        inputs = []
        for index in range(len(genes)):
            inputs.append(dict(zip(self.symbols[0][1], terminalInputs)))

        trees = []
        for index in range(len(genes)):
            gene = genes[index][0]
            if isinstance(gene[0], int):
                node = Node(gene[0], 0)
            else:
                if genes[index][1]:
                    node = Node(gene[0], self.genomes[1].arities[gene[0]])
                else:
                    node = Node(gene[0], self.genomes[0].arities[gene[0]])

            trees.append((node, genes[index][1]))
            for symbol in gene[1:]:
                if isinstance(symbol, int):
                    node.append(symbol, 0)
                else:
                    if genes[index][1]:
                        node.append(symbol, self.genomes[1].arities[symbol])
                    else:
                        node.append(symbol, self.genomes[0].arities[symbol])

        homeoticInputs = {}
        for index in range(len(trees)):
            if trees[index][1]:
                homeoticInputs[index] = trees[index][0].eval(self.genomes[1], homeoticInputs)
            else:
                homeoticInputs[index] = trees[index][0].eval(self.genomes[0], inputs[index])

        return homeoticInputs

    def copy(self):
        new = Chromosome(None, None, None, False)
        new.fitness = self.fitness
        new.genomes = self.genomes
        new.head = self.head
        new.headRanges = self.headRanges
        new.tailRanges = self.tailRanges
        new.homeoticRanges = self.homeoticRanges
        new.sequence = self.sequence[:]
        new.symbols = self.symbols
        return new

    def mutation(self, rate):
        indices = []
        for index in range(len(self.headRanges)):
            indices += list(zip(list(range(self.headRanges[index][0], self.tailRanges[index][1])), [index for i in range(self.headRanges[index][0], self.tailRanges[index][1])]))

        random.shuffle(indices)
        indices = indices[:int(round(rate * len(indices)))]
        for index in indices:
            if self.head[index[0]]:
                if self.homeoticRanges[index[0]]:
                    symbols = self.symbols[1][0] + self.symbols[1][1]
                else:
                    symbols = self.symbols[0][0] + self.symbols[0][1]

                symbols.remove(self.sequence[index[0]])
                if len(symbols) > 0:
                    self.sequence[index[0]] = random.choice(symbols)

            else:
                if self.homeoticRanges[index[0]]:
                    symbols = self.symbols[1][1][:]
                else:
                    symbols = self.symbols[0][1][:]

                symbols.remove(self.sequence[index[0]])
                if len(symbols) > 0:
                    self.sequence[index[0]] = random.choice(symbols)

    def inversion(self):
        index = random.randrange(len(self.headRanges))
        start = random.randrange(self.headRanges[index][0], self.headRanges[index][1])
        end = random.randint(start + 1, self.headRanges[index][1])
        self.sequence[start:end] = reversed(self.sequence[start:end])

    def ISTransposition(self):
        index = random.randrange(len(self.headRanges))
        start = random.randrange(self.headRanges[index][0], self.tailRanges[index][1])
        end = random.randint(start + 1, self.tailRanges[index][1])
        sequence = self.sequence[start:end]
        destIndex = 0
        while True:
            destIndex = random.randrange(len(self.headRanges))
            if self.homeoticRanges[index] is self.homeoticRanges[destIndex]:
                break

        destination = random.randrange(self.headRanges[destIndex][0] + 1, self.headRanges[destIndex][1])
        destSequence = self.sequence[destination:self.headRanges[destIndex][1]]
        sequence += destSequence
        self.sequence[destination:self.headRanges[destIndex][1]] = sequence[:self.headRanges[destIndex][1] - destination]

    def RISTransposition(self, tries=0):
        if tries > sys.getrecursionlimit() - 5:
            return

        index = random.randrange(len(self.headRanges))
        start = random.randrange(self.headRanges[index][0], self.headRanges[index][1])
        for i in range(start, self.headRanges[index][1]):
            if self.homeoticRanges[index]:
                if self.sequence[i] in self.symbols[1][0]:
                    start = i
                    if self.sequence[i] not in self.symbols[1][0]:
                        self.RISTransposition(tries + 1)
                        return

            else:
                if self.sequence[i] in self.symbols[0][0]:
                    start = i
                    if self.sequence[i] not in self.symbols[0][0]:
                        self.RISTransposition(tries + 1)
                        return

        end = random.randint(start + 1, self.tailRanges[index][1])
        sequence = self.sequence[start:end]
        destIndex = 0
        while True:
            destIndex = random.randrange(len(self.headRanges))
            if self.homeoticRanges[index] is self.homeoticRanges[destIndex]:
                break

        destination = self.headRanges[destIndex][0]
        destSequence = self.sequence[destination:self.headRanges[destIndex][1]]
        sequence += destSequence
        self.sequence[destination:self.headRanges[destIndex][1]] = sequence[:self.headRanges[destIndex][1] - destination]

    def geneTransposition(self):
        index = random.randrange(len(self.headRanges))
        destIndex = 0
        while True:
            destIndex = random.randrange(len(self.headRanges))
            if self.homeoticRanges[index] is self.homeoticRanges[destIndex]:
                if self.headRanges[index][1] - self.headRanges[index][0] == self.headRanges[destIndex][1] - self.headRanges[destIndex][0]:
                    break

        self.sequence[self.headRanges[destIndex][0]:self.tailRanges[destIndex][1]] = self.sequence[self.headRanges[index][0]:self.tailRanges[index][1]]

    def onePointRecombination(self, other):
        point = random.randrange(len(self.sequence))
        self.sequence[point:], other.sequence[point:] = other.sequence[point:], self.sequence[point:]

    def twoPointRecombination(self, other):
        firstPoint = random.randrange(len(self.sequence))
        secondPoint = random.randrange(len(self.sequence))
        self.sequence[firstPoint:secondPoint], other.sequence[firstPoint:secondPoint] = other.sequence[firstPoint:secondPoint], self.sequence[firstPoint:secondPoint]

    def geneRecombination(self, other):
        index = random.randrange(len(self.headRanges))
        firstPoint = self.headRanges[index][0]
        secondPoint = self.tailRanges[index][1]
        self.sequence[firstPoint:secondPoint], other.sequence[firstPoint:secondPoint] = other.sequence[firstPoint:secondPoint], self.sequence[firstPoint:secondPoint]

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
        self.arities[symbol] = 0

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

    def symbols(self, numGenes, homeotic):
        if homeotic:
            return (list(self.functions.keys()), list(range(numGenes)))

        return (list(self.functions.keys()), self.terminals)


    def _add(inputs):
        if inputs[0] == "Error" or inputs[1] == "Error":
            return "Error"
        return inputs[0] + inputs[1]
    def _sub(inputs):
        if inputs[0] == "Error" or inputs[1] == "Error":
            return "Error"
        return inputs[0] - inputs[1]
    def _mul(inputs):
        if inputs[0] == "Error" or inputs[1] == "Error":
            return "Error"
        return inputs[0] * inputs[1]
    def _div(inputs):
        if inputs[0] == "Error" or inputs[1] == "Error":
            return "Error"
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
    def __init__(self, symbol, arity):
        self.children = []
        self.symbol = symbol
        self.arity = arity

    def append(self, symbol, arity):
        if len(self.children) < self.arity:
            self.children.append(Node(symbol, arity))
            return True
        else:
            for child in self.children:
                if child.append(symbol, arity):
                    return True

            return False

    def eval(self, genome, terminalInputs):
        inputs = []
        for child in self.children:
            inputs.append(child.eval(genome, terminalInputs))

        if self.symbol in list(genome.functions.keys()):
            return genome.functions[self.symbol](inputs)

        return terminalInputs[self.symbol]


genome = Genome()
genome.addFunction("A", Genome._AND, 2)
genome.addFunction("O", Genome._OR, 2)
genome.addFunction("N", Genome._NOT, 1)
genome.addTerminal("a")
genome.addTerminal("b")
chromosome = Chromosome([2, 3, 4], [], [genome, genome])
other = Chromosome([2, 3, 4], [], [genome, genome])
print(chromosome.headRanges)
print(chromosome.tailRanges)
print(chromosome.homeoticRanges)
print(chromosome.symbols)
print(chromosome.head)
print(chromosome.sequence)
print(chromosome.eval([True, False]))