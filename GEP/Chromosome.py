__author__ = 'Coleman'

from GEP.Gene import Gene
import random

class Chromosome():
    def __init__(self):
        self.genes = []
        self.homeotics = []

    def initRand(self, numGenes, numHomeotics, headLength, genome):
        self.genes = [Gene(genome, False).initRand(headLength, numGenes) for i in range(numGenes)]
        self.homeotics = [Gene(genome, True).initRand(headLength, numGenes) for i in range(numHomeotics)]

    def replicate(self):
        newGenes = []
        newHomeotics = []
        for gene in self.genes:
            newGenes.append(gene.replicate())

        for homeotic in self.homeotics:
            newHomeotics.append(homeotic.replicate())

        newChromosome = Chromosome()
        newChromosome.genes = newGenes
        newChromosome.homeotics = newHomeotics

        return newChromosome

    def eval(self, inputs):
        returnValues = []
        for gene in self.genes:
            evalInputs = {terminal: value for (terminal, value) in list(zip(gene.genome.terminals, inputs))}
            returnValues.append(gene.eval(evalInputs))

        if len(self.homeotics) > 0:
            homeoticReturnValues = []
            outputs = {index: output for (index, output) in list(zip(range(0, len(self.genes)), returnValues))}
            for homeotic in self.homeotics:
                homeoticReturnValues.append(homeotic.eval(outputs))

            returnValues = homeoticReturnValues

        return returnValues

    def mutation(self, rate, homeoticRate):
        for gene in self.genes:
            head = gene.head
            tail = gene.tail
            functions = list(gene.genome.functions.keys())
            functions.pop(functions.index("arity"))
            for i in range(len(head)):
                if random.random() < rate:
                    head[i] = random.choice(functions + gene.genome.terminals)

            for i in range(len(tail)):
                if random.random() < rate:
                    tail[i] = random.choice(gene.genome.terminals)

        rate *= homeoticRate
        for homeotic in self.homeotics:
            head = homeotic.head
            tail = homeotic.tail
            functions = list(homeotic.genome.functions.keys())
            functions.pop(functions.index("arity"))
            for i in range(len(head)):
                if random.random() < rate:
                    head[i] = random.choice(functions + list(range(len(self.genes))))

            for i in range(len(tail)):
                if random.random() < rate:
                    tail[i] = random.randint(0, len(self.genes) - 1)

    def inversion(self, rate, homeoticRate):
        for gene in self.genes:
            if random.random() < rate:
                start = random.randint(0, len(gene.head) - 1)
                stop = random.randint(start, len(gene.head) - 1)
                gene.head[start:stop + 1] = reversed(gene.head[start:stop + 1])

        rate *= homeoticRate
        for homeotic in self.homeotics:
            if random.random() < rate:
                start = random.randint(0, len(homeotic.head) - 1)
                stop = random.randint(start, len(homeotic.head) - 1)
                homeotic.head[start:stop + 1] = reversed(homeotic.head[start:stop + 1])

    def ISTransposition(self, rate, homeoticRate):
        if random.random() < rate:
            sequence = []
            for gene in self.genes:
                sequence = sequence + gene.head + gene.tail

            gene = random.choice(self.genes)
            headLength = len(gene.head)
            start = random.randint(0, len(sequence) - 1)
            stop = random.randint(start, start + headLength)
            ISSequence = sequence[start:stop + 1]
            insertionPoint = random.randint(1, headLength - 1)
            for i in range(len(ISSequence)):
                gene.head.insert(insertionPoint + i, ISSequence[i])

            gene.head = gene.head[:headLength]

        rate *= homeoticRate
        if random.random() < rate:
            sequence = []
            for homeotic in self.homeotics:
                sequence = sequence + homeotic.head + homeotic.tail

            homeotic = random.choice(self.homeotics)
            headLength = len(homeotic.head)
            start = random.randint(0, len(sequence) - 1)
            stop = random.randint(start, start + headLength)
            ISSequence = sequence[start:stop + 1]
            insertionPoint = random.randint(1, headLength - 1)
            for i in range(len(ISSequence)):
                homeotic.head.insert(insertionPoint + i, ISSequence[i])

            homeotic.head = homeotic.head[:headLength]

    def RISTransposition(self, rate, homeoticRate):
        if random.random() < rate:
            sequence = []
            for gene in self.genes:
                sequence = sequence + gene.head + gene.tail

            gene = random.choice(self.genes)
            headLength = len(gene.head)
            functions = list(gene.genome.functions.keys())
            functions.pop(functions.index("arity"))
            start = random.randint(0, len(sequence) - 1)
            while sequence[start] not in functions:
                start += 1
                if start == len(sequence) - 1:
                    break

            stop = random.randint(start, start + headLength)
            RISSequence = sequence[start:stop + 1]
            for i in range(len(RISSequence)):
                gene.head.insert(i, RISSequence[i])

            gene.head = gene.head[:headLength]

        rate *= homeoticRate
        if random.random() < rate:
            sequence = []
            for homeotic in self.homeotics:
                sequence = sequence + homeotic.head + homeotic.tail

            homeotic = random.choice(self.homeotics)
            headLength = len(homeotic.head)
            functions = list(gene.genome.functions.keys())
            functions.pop(functions.index("arity"))
            start = random.randint(0, len(sequence) - 1)
            while sequence[start] not in functions:
                start += 1
                if start == len(sequence) - 1:
                    break

            stop = random.randint(start, start + headLength)
            RISSequence = sequence[start:stop + 1]
            for i in range(len(RISSequence)):
                homeotic.head.insert(i, RISSequence[i])

            homeotic.head = homeotic.head[:headLength]

    def geneTransposition(self, rate, homeoticRate):
        if random.random() < rate:
            self.genes[0] = random.choice(self.genes).replicate()

        rate *= homeoticRate
        if random.random() < rate:
            self.homeotics[0] = random.choice(self.homeotics).replicate()

    def onePointRecombination(self, rate, homeoticRate, otherChromosome):
        if random.random() < rate:
            sequence = []
            otherSequence = []
            for gene in self.genes:
                sequence = sequence + gene.head + gene.tail

            for otherGene in otherChromosome.genes:
                otherSequence = otherSequence + otherGene.head + otherGene.tail

            if len(sequence) != len(otherSequence):
                return

            recombinationPoint = random.randint(0, len(sequence) - 1)
            sequence[:recombinationPoint], otherSequence[recombinationPoint:] = otherSequence[:recombinationPoint], sequence[recombinationPoint:]
            for gene in self.genes:
                gene.head = sequence[:len(gene.head)]
                for i in range(len(gene.head)):
                    sequence.pop(0)

                gene.tail = sequence[:len(gene.tail)]
                for i in range(len(gene.tail)):
                    sequence.pop(0)

            for otherGene in otherChromosome.genes:
                otherGene.head = otherSequence[:len(otherGene.head)]
                for i in range(len(otherGene.head)):
                    otherSequence.pop(0)

                otherGene.tail = sequence[:len(otherGene.tail)]
                for i in range(len(otherGene.tail)):
                    otherSequence.pop(0)

        rate *= homeoticRate
        if random.random() < rate:
            sequence = []
            otherSequence = []
            for homeotic in self.homeotics:
                sequence = sequence + homeotic.head + homeotic.tail

            for otherHomeotic in otherChromosome.homeotics:
                otherSequence = otherSequence + otherHomeotic.head + otherHomeotic.tail

            if len(sequence) != len(otherSequence):
                return

            recombinationPoint = random.randint(0, len(sequence) - 1)
            sequence[:recombinationPoint], otherSequence[recombinationPoint:] = otherSequence[:recombinationPoint], sequence[recombinationPoint:]
            for homeotic in self.homeotics:
                homeotic.head = sequence[:len(homeotic.head)]
                for i in range(len(homeotic.head)):
                    sequence.pop(0)

                homeotic.tail = sequence[:len(homeotic.tail)]
                for i in range(len(homeotic.tail)):
                    sequence.pop(0)

            for otherHomeotic in otherChromosome.homeotics:
                otherHomeotic.head = otherSequence[:len(otherHomeotic.head)]
                for i in range(len(otherHomeotic.head)):
                    otherSequence.pop(0)

                otherHomeotic.tail = sequence[:len(otherHomeotic.tail)]
                for i in range(len(otherHomeotic.tail)):
                    otherSequence.pop(0)

    def twoPointRecombination(self, rate, homeoticRate, otherChromosome):
        if random.random() < rate:
            sequence = []
            otherSequence = []
            for gene in self.genes:
                sequence = sequence + gene.head + gene.tail

            for otherGene in otherChromosome.genes:
                otherSequence = otherSequence + otherGene.head + otherGene.tail

            if len(sequence) != len(otherSequence):
                return

            recombinationPoint = random.randint(0, len(sequence) - 1)
            otherPoint = random.randint(recombinationPoint, len(sequence) - 1)
            sequence[:recombinationPoint], otherSequence[recombinationPoint:otherPoint], sequence[otherPoint:] = otherSequence[:recombinationPoint], sequence[recombinationPoint:otherPoint], sequence[otherPoint:]
            for gene in self.genes:
                gene.head = sequence[:len(gene.head)]
                for i in range(len(gene.head)):
                    sequence.pop(0)

                gene.tail = sequence[:len(gene.tail)]
                for i in range(len(gene.tail)):
                    sequence.pop(0)

            for otherGene in otherChromosome.genes:
                otherGene.head = otherSequence[:len(otherGene.head)]
                for i in range(len(otherGene.head)):
                    otherSequence.pop(0)

                otherGene.tail = sequence[:len(otherGene.tail)]
                for i in range(len(otherGene.tail)):
                    otherSequence.pop(0)

        rate *= homeoticRate
        if random.random() < rate:
            sequence = []
            otherSequence = []
            for homeotic in self.homeotics:
                sequence = sequence + homeotic.head + homeotic.tail

            for otherHomeotic in otherChromosome.homeotics:
                otherSequence = otherSequence + otherHomeotic.head + otherHomeotic.tail

            if len(sequence) != len(otherSequence):
                return

            recombinationPoint = random.randint(0, len(sequence) - 1)
            otherPoint = random.randint(recombinationPoint, len(sequence) - 1)
            sequence[:recombinationPoint], otherSequence[recombinationPoint:] = otherSequence[:recombinationPoint], sequence[recombinationPoint:]
            for homeotic in self.homeotics:
                homeotic.head = sequence[:len(homeotic.head)]
                for i in range(len(homeotic.head)):
                    sequence.pop(0)

                homeotic.tail = sequence[:len(homeotic.tail)]
                for i in range(len(homeotic.tail)):
                    sequence.pop(0)

            for otherHomeotic in otherChromosome.homeotics:
                otherHomeotic.head = otherSequence[:len(otherHomeotic.head)]
                for i in range(len(otherHomeotic.head)):
                    otherSequence.pop(0)

                otherHomeotic.tail = sequence[:len(otherHomeotic.tail)]
                for i in range(len(otherHomeotic.tail)):
                    otherSequence.pop(0)

    def geneRecombination(self, rate, homeoticRate, otherGene):
        if random.random() < rate:
            if len(self.genes) != len(otherChromosome.genes):
                return
            
            recombinationPoint = random.randint(0, len(self.genes) - 1)
            self.genes[:recombinationPoint], otherChromosome.genes[recombinationPoint:] = otherChromosome.genes[:recombinationPoint], self.genes[recombinationPoint:]
            
        rate *= homeoticRate
        if random.random() < rate:
            if len(self.homeotics) != len(otherChromosome.homeotics):
                return
            
            recombinationPoint = random.randint(0, len(self.homeotics) - 1)
            self.homeotics[:recombinationPoint], otherChromosome.homeotics[recombinationPoint:] = otherChromosome.homeotics[:recombinationPoint], self.homeotics[recombinationPoint:]
            
    def printChromosome(self):
        for gene in self.genes:
            print(gene.head + gene.tail, end=':\n')

        for homeotic in self.homeotics:
            print(homeotic.head + homeotic.tail, end=';\n')

        print()
