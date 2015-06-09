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

class Chromosome:
    def __init__(self, headLengths, homeoticLengths, genomes, gen=True):
        self.fitness = 0
        self.terminals = None
        self.symbols = None
        self.genomes = None

        self.headRanges = None
        self.tailRanges = None
        self.homeoticRanges = None
        self.sequence = None

        self.head = None
        if gen:
            terminals = [genomes[i].terminals for i in range(len(headLengths))]
            homeoticTerminals = [list(range(len(headLengths))) for i in range(len(homeoticLengths))]
            self.terminals = terminals + homeoticTerminals
            self.symbols = [genomes[i].symbols(len(headLengths), i) for i in range(len(genomes))]
            self.genomes = genomes

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
            inputs.append(dict(zip(self.symbols[index][1], terminalInputs)))

        trees = []
        for index in range(len(genes)):
            gene = genes[index][0]
            if isinstance(gene[0], int):
                node = Node(gene[0], 0)
            else:
                node = Node(gene[0], self.genomes[index].arities[gene[0]])

            trees.append((node, genes[index][1]))
            for symbol in gene[1:]:
                if isinstance(symbol, int):
                    node.append(symbol, 0)
                else:
                    node.append(symbol, self.genomes[index].arities[symbol])

        homeoticInputs = {}
        for index in range(len(trees)):
            if trees[index][1]:
                homeoticInputs[index] = trees[index][0].eval(self.genomes[index], homeoticInputs)
            else:
                homeoticInputs[index] = trees[index][0].eval(self.genomes[index], inputs[index])

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
        new.terminals = self.terminals
        return new

    def mutation(self, rate):
        indices = []
        for index in range(len(self.headRanges)):
            indices += list(zip(list(range(self.headRanges[index][0], self.tailRanges[index][1])), [index for i in range(self.headRanges[index][0], self.tailRanges[index][1])]))

        random.shuffle(indices)
        indices = indices[:int(round(rate * len(indices)))]
        for index in indices:
            if self.head[index[0]]:
                symbols = self.symbols[index[1]][0] + self.symbols[index[1]][1]
                symbols.remove(self.sequence[index[0]])
                if len(symbols) > 0:
                    self.sequence[index[0]] = random.choice(symbols)
            else:
                symbols = self.symbols[index[1]][1][:]
                symbols.remove(self.sequence[index[0]])
                if len(symbols) > 0:
                    self.sequence[index[0]] = random.choice(symbols)

    def inversion(self):
        index = random.randint(0, len(self.headRanges) - 1)
        start = random.randint(self.headRanges[index][0], self.headRanges[index][1] - 1)
        end = random.randint(start + 1, self.headRanges[index][1])
        print(start, end)
        self.sequence[start:end] = reversed(self.sequence[start:end])

    def ISTransposition(self):
        index = random.randint(0, len(self.headRanges) - 1)
        start = random.randint(self.headRanges[index][0], self.tailRanges[index][1] - 1)
        end = random.randint(start + 1, self.tailRanges[index][1])
        sequence = self.sequence[start:end]
        destIndex = 0
        while True:
            destIndex = random.randint(0, len(self.headRanges) - 1)
            if self.homeoticRanges[index] is self.homeoticRanges[destIndex]:
                break

        destination = random.randint(self.headRanges[destIndex][0] + 1, self.headRanges[destIndex][1] - 1)
        if self.headRanges[destIndex][1] - destination < end - start:
            self.sequence[destination:self.headRanges[destIndex][1]] = sequence[:self.headRanges[destIndex][1] - destination]
        else:
            self.sequence[destination:destination + (end - start)] = sequence

        print(start, end)
        print(destination)

    def RISTransposition(self):
        pass

    def geneTransposition(self):
        pass

    def onePointRecombination(self, other):
        pass

    def twoPointRecombination(self, other):
        pass

    def geneRecombination(self, other):
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

    def symbols(self, numGenes, index):
        if index > numGenes - 1:
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
chromosome = Chromosome([5, 5, 5], [5], [genome for i in range(4)])
print(chromosome.headRanges)
print(chromosome.tailRanges)
print(chromosome.homeoticRanges)
print(chromosome.symbols)
print(chromosome.head)
print(chromosome.sequence)
chromosome.ISTransposition()
print(chromosome.sequence)