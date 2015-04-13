__author__ = 'Coleman'

from Gene import Gene

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
        
    def mutation():
        pass
    
    def inversion():
        pass
    
    def ISTransposition():
        pass
    
    def RISTransposition():
        pass
    
    def geneTransposition():
        pass
    
    def onePointRecombination():
        pass
    
    def twoPointRecombination():
        pass
    
    def geneRecombination():
        pass
