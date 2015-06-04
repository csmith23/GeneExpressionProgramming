__author__ = 'Coleman'

from GEP.Chromosome import Chromosome
import random

class Environment:
    def __init__(self):
        self.population = []

        self.homeoticRate = 0 # multiplied by the other rates for every homeotic gene
        self.mutationRate = 0
        self.inversionRate = 0
        self.ISTranspositionRate = 0
        self.RISTranspositionRate = 0
        self.geneTranspositionRate = 0
        self.onePointRecombinationRate = 0
        self.twoPointRecombinationRate = 0
        self.geneTranspositionRate = 0

    def init(self, populationSize, numGenes, numHomeotics, headLength, homeoticHeadLength, genome):
        self.population = [Chromosome().initRand(numGenes, numHomeotics, headLength, homeoticHeadLength, genome) for i in range(populationSize)]

    def setRates(self, homeoticRate=0, mutationRate=0, inversionRate=0, ISTranspositionRate=0, RISTranspositionRate=0, geneTranspositionRate=0, onePointRecombinationRate=0, twoPointRecombinationRate=0, geneRecombinationRate=0):
        self.homeoticRate = homeoticRate
        self.mutationRate = mutationRate
        self.inversionRate = inversionRate
        self.ISTranspositionRate = ISTranspositionRate
        self.RISTranspositionRate = RISTranspositionRate
        self.geneTranspositionRate = geneTranspositionRate
        self.onePointRecombinationRate = onePointRecombinationRate
        self.twoPointRecombinationRate = twoPointRecombinationRate
        self.geneRecombinationRate = geneRecombinationRate

    def run(self, inputsOutputs, fitnessFunction):
        '''fitness must have a argument list of ([actual], [target], max=False)
        inputsOutputs : (inputs, outputs)
        actual : []*
        target : []*
        if max = True, fitness must return the maximum fitness value of the function
        fitness must return an integer value, so if it would return a percentage, make it an integer with a maximum of 100
        '''
        self.printChromosomes(0)
        generation = 0
        while True:
            generation += 1
            for i in range(len(self.population)):
                chromosome = self.population[i]
                answers = []
                for input in inputsOutputs[0]:
                    answers.append(chromosome.eval(input))

                fitness = fitnessFunction(answers, inputsOutputs[1])
                if fitness == fitnessFunction([], [], True):
                    chromosome.printChromosome()
                    return chromosome

                self.population[i] = (chromosome, fitness)

            selected = self.select()
            self.replicate(selected)
            self.modify()
            self.printChromosomes(generation)

    def select(self):
        chromosomes = []
        bestFitness = 0
        bestIndex = 0
        for i in range(len(self.population)):
            fitness = self.population[i][1]
            if fitness > bestFitness:
                bestIndex = i

            for j in range(fitness):
                chromosomes.append(self.population[i][0])

        selected = []
        for i in range(len(self.population) - 1):
            index = random.randint(0, len(chromosomes) - 1)
            selected.append(chromosomes[index])

        selected.append(self.population[bestIndex][0])

        return selected

    def replicate(self, chromosomes):
        self.population = []
        for chromosome in chromosomes:
            self.population.append(chromosome.replicate())

    def modify(self):
        for i in range(len(self.population)):
            chromosome = self.population[i]
            chromosome.mutation(self.mutationRate, self.homeoticRate)
            chromosome.inversion(self.inversionRate, self.homeoticRate)
            chromosome.ISTransposition(self.ISTranspositionRate, self.homeoticRate)
            chromosome.RISTransposition(self.RISTranspositionRate, self.homeoticRate)
            chromosome.geneTransposition(self.geneTranspositionRate, self.homeoticRate)
            otherIndex = i
            while otherIndex == i:
                otherIndex = random.randint(0, len(self.population) - 1)

            chromosome.onePointRecombination(self.onePointRecombinationRate, self.homeoticRate, self.population[otherIndex])
            otherIndex = i
            while otherIndex == i:
                otherIndex = random.randint(0, len(self.population) - 1)

            chromosome.twoPointRecombination(self.twoPointRecombinationRate, self.homeoticRate, self.population[otherIndex])
            otherIndex = i
            while otherIndex == i:
                otherIndex = random.randint(0, len(self.population) - 1)

            chromosome.geneRecombination(self.geneRecombinationRate, self.homeoticRate, self.population[otherIndex])

    def printChromosomes(self, generation):
        print(generation)
        for i in range(len(self.population)):
            print(i, ':', end='')
            for gene in self.population[i].genes:
                print(gene.head + gene.tail, end=':')

            for homeotic in self.population[i].homeotics:
                print(homeotic.head + homeotic.tail, end=':')

            print()

        print()