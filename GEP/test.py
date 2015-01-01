__author__ = 'Coleman'

from GEP.Genome import Genome
from GEP.Chromosome import *

# genome = Genome()
#
# genome.functions.update(Genome.ARITHMETIC_SET)
# genome.terminals = ["a", "b", "c"]
#
# genes = []
# for i in range(3):
#     newGene = Gene(genome)
#     newGene.initRand(7)
#     genes.append(newGene)
#
# link = Gene(genome)
# link.head = ["+", "+", "0", "1", "2", "2", "1"]
# link.tail = ["2", "2", "1", "0", "1", "2", "2", "1"]
# genes.append(link)
# chromosome = Chromosome(genome, 1, 1)
# chromosome.genes = genes
# tree = chromosome.Tree()
# chromosome.print()
# termdict = {"a": 12, "b": 3, "c": 6}
# output = tree.step(termdict)
# print(output)

# NODE----------------------------------------
    # def print(self, depth = 0):
    #     '''a method of printing a whole tree recursively, post order
    #     -TERM means the preceding output came from a leaf node'''
    #     if len(self.children) > 0:
    #         # for child in self.children:
    #         #     child.print(depth + 1)
    #         # self.children[0].print(depth + 1)
    #         print(depth, end = ":")
    #         print(self.character)
    #         # self.children[0].print(depth + 1)
    #         for child in self.children:
    #             child.print(depth + 1)
    #     else:
    #         print(depth, end = ":")
    #         print(self.character, end = ":")
    #         print("-TERM")
# GENE-----------------------------------------
    # def print(self):
    #     print(self.head, end = ":")
    #     print(self.tail)
    #     print()
# CHROMOSOME------------------------------------
    # def print(self):
    #     for gene in self.genes:
    #         gene.print()
    #
    #     tree = self.Tree()
    #     tree[-1].print()


# 2 genetic operators, p91

# Replication, exact copying
# Selection, roulette wheel, with optional elitism

# 8 modification operators, p91; execute on each chromosome in this order

# within one gene in a chromosome-
# Mutation, head and tail; randomly alters a single element in the head or tail
#     random new character, run for every character in the chromosome
#     return 1
# Inversion, head; inverts a sequence of random length in the head
#     random starting point and length
#     return 1

# within one chromosome, between or not between genes-
# IS transposition, head; any sequence to anywhere in the head except the root
#     random index, length, other gene, and new index
#     return 1
# RIS transposition, head; any sequence beginning with a function to anywhere in the head: uses scanning from index, with nop on fail
#     random index, length, other gene, and new index
#     return 1
# Gene transposition, head and tail; any gene overriding any gene, even itself
#     random gene, other gene
#     return 1

# within one setting, between chromosomes-
# One-point recombination, head
#     random index
#     return 2
# Two-point recombination, head
#     random index, length
#     return 2
# Gene recombination, head and tail
#     random gene index
#     return 2

# use each gene concatenated and slice notation to return them as separate, even in head and tail