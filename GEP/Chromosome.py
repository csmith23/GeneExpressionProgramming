__author__ = 'Coleman'

import random
import copy
from collections import deque

class Chromosome():
    def __init__(self, genome):
        self.genome = genome
        self.fitness = 0
        self.genes = []
        self.link = Gene()
        self.tree = Node()
        self.headLength = 0
        # self.genes.append(Gene(genome).initRand(length) for i in range(numGenes))

    def initRand(self, numGenes, length):
        self.genes.append(Gene(self.genome).initRand(length) for i in range(numGenes))
        self.headLength = length

    def Tree(self):
        subTrees = []
        for gene in self.genes:
            subTrees.append(gene.getTree(subTrees))

        self.tree = self.link.getTree(subTrees)

    def getSequence(self):
        sequence = []
        for gene in self.genes:
            sequence.append(gene.getSequence())

        return sequence

    def mutation(self, rate):
        symbols = self.genome.terminals + list(self.genome.functions.keys())
        symbols.remove("arity")

        for gene in self.genes:
            for i in range(len(gene.head)):
                if random.random() < rate:
                    gene.head[i] = random.choice(symbols)
            for i in range(len(gene.tail)):
                if random.random < rate:
                    gene.tail[i] = random.choice(self.genome.terminals)

    def inversion(self, rate):
        if random.random() < rate:
            gene = random.choice(self.genes)
            index = random.randint(len(gene.head))
            length = random.randint(1, len(gene.head) - index)
            gene.head[index:index + length] = reversed(gene.head[index:index + length])

    def ISTransposition(self, rate):
        '''pick index, length, then new gene index and index not 0, then truncate with this = this[:headLength]'''
        if random.random() < rate:
            sequence = self.getSequence()
            index = random.randint(0, len(sequence))
            length = random.randint(1, len(sequence) - index)
            destinationGene = random.choice(self.genes)
            destinationIndex = random.randint(1, self.headLength)
            IS = sequence[index:length]
            destinationGene.head.insert(destinationIndex, IS)
            del destinationGene.head[self.headLength:]

    def RISTransposition(self, rate):
        '''pick index, length, then new gene index and index, find first function, then truncate in same manner'''
        if random.random() < rate:
            sequence = self.getSequence()
            index = random.randint(0, len(sequence))
            while sequence[index] not in self.genome.keys():
                index += 1
                if index == len(sequence):
                    return

            length = random.randint(1, len(sequence) - index)
            destinationGene = random.choice(self.genes)
            IS = sequence[index:length]
            destinationGene.head.insert(0, IS)
            del destinationGene.head[self.headLength:]

    def geneTransposition(self, rate):
        '''pick gene index, other gene index'''
        if random.random() < rate:
            index = random.choice(self.genes)
            while True:
                destinationIndex = random.randint(0, len(self.genes))
                if index != destinationIndex:
                    break

            self.genes[destinationIndex] = copy.deepcopy(self.genes[index])

    '''all recombination functions must first turn the gene list into one string, then back into separate genes'''
    def onePointRecombination(self, rate, other):
        '''choose other chromosome(should be the chromosome at i + 1 or -1 if i == length - 1), index in range 1 to length of chromosome - 1, then recombine with this[:index] + other[index:] and vice versa'''
        if random.random() < rate:
            sequence = self.getSequence()
            otherSequence = other.getSequence()
            if len(sequence) != len(otherSequence):
                return

            index = random.randint(1, len(sequence) - 1)

            sequence, otherSequence = sequence[:index] + otherSequence[index:], otherSequence[:index] + sequence[index:]

            for chromosome, sequence in [(self, sequence), (other, otherSequence)]:
                for gene in chromosome.genes:
                    gene.head, sequence = sequence[:len(gene.head)], sequence[len(gene.head):]
                    gene.tail, sequence = sequence[:len(gene.tail)], sequence[len(gene.tail):]

    def twoPointRecombination(self, rate, other):
        '''choose other chromosome(should be the chromosome at i + 1 or -1 if i == length - 1), two indices in range 1 to length of chromosome - 1, then recombine with this[:index1] + other[index1:index2] + this[index2:] and vice versa'''
        if random.random() < rate:
            sequence = self.getSequence()
            otherSequence = other.getSequence()
            if len(sequence) != len(otherSequence):
                return

            index = random.randint(1, len(sequence) - 1)
            while True:
                otherIndex = random.randint(1, len(sequence) - 1)
                if otherIndex != index:
                    break
            index, otherIndex = random.sort([index, otherIndex])

            sequence, otherSequence = sequence[:index] + otherSequence[index:otherIndex] + sequence[otherIndex:], otherSequence[:index] + sequence[index:otherIndex] + otherSequence[otherIndex:]

            for chromosome, sequence in [(self, sequence), (other, otherSequence)]:
                for gene in chromosome.genes:
                    gene.head, sequence = sequence[:len(gene.head)], sequence[len(gene.head):]
                    gene.tail, sequence = sequence[:len(gene.tail)], sequence[len(gene.tail):]

    def geneRecombination(self, rate, other):
        '''choose other chromosome(should be the chromosome at i + 1 or -1 if i == length - 1), a gene index, and swap the two genes at that index'''
        if random.random() < rate:
            if len(self.genes) != len(other.genes):
                return

            index = random.randint(0, len(self.genes) - 1)
            self.genes, other.genes = self.genes[:index] + other.genes[index:index + 1] + self.genes[index + 1:], other.genes[:index] + self.genes[index:index + 1] + other.genes[index + 1:]

class Gene():
    def __init__(self, genome):
        self.genome = genome
        self.head = []
        self.tail = []

    def initRand(self, length):
        symbols = self.genome.terminals + list(self.genome.functions.keys())
        symbols.remove("arity")

        for i in range(length):
            self.head.append(random.choice(symbols))

        for i in range(length * (self.genome.functions["arity"] - 1) + 1):
            self.tail.append(random.choice(self.genome.terminals))

    def getSequence(self):
        return self.head + self.tail

    def getTree(self, trees):
        sequence = deque(self.getSequence())
        queue = deque()
        root = Node(sequence.popleft(), self.genome)
        current = root
        for character in sequence:
            print(character)
            while current.character in self.genome.terminals or len(current.children) == self.genome.functions[current.character][1]:
                if len(queue) == 0:
                    return root
                current = queue.popleft()
            if character.isdigit():
                current.children.append(trees[int(float(character))])
            else:
                child = Node(character, self.genome)
                current.children.append(child)
                queue.append(child)

        return root

class Node():
    def __init__(self, character, genome):
        self.genome = genome
        self.parent = None
        self.children = []
        self.character = character

    def step(self, terminalDict):
        inputs = []
        if self.character in list(self.genome.functions.keys()):
            for child in self.children:
                inputs.append(child.step(terminalDict))
            return self.genome.functions[self.character][0](inputs)
        else:
            return terminalDict[self.character]
