__author__ = 'Coleman'

from GEP.Chromosome import Chromosome

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
        self.functions = {} #dictionary - character : (function, arity)
        self.terminals = [] #list - character
        self.arity = 0

    def addFunctions(self, otherDict):
        otherDict = dict(otherDict)
        if ("arity" in self.functions) and ("arity" in otherDict):
            if otherDict["arity"] < self.functions["arity"]:
                otherDict["arity"] = self.functions["arity"]
        self.function.update(dict)

    def Chromosome(self, numGenes,  length):
        return Chromosome(self, numGenes, length)