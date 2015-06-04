__author__ = 'Coleman'

import random

class Gene():
    def __init__(self, genome):
        self.genome = genome
        self.head = []
        self.tail = []
        self.elementLayers = []

    def initRand(self, length):
        arity = self.genome.functions["arity"]
        functions = self.genome.functions.keys()
        functions.remove("arity")
        terminals = self.genome.terminals
        self.head = [random.choice(functions + terminals) for i in range(length)]
        self.tail = [random.choice(terminals) for i in range(length * (arity - 1) + 1)]
        return self

    def replicate(self):
        newGene = Gene(self.genome)
        newGene.head = self.head[:]
        self.tail = self.tail[:]
        return newGene

    def eval(self, inputs, chromosome):
        self.elementLayers = []

        self.orderStack()
        evalStack = self.evalRecur()
        return self.evalStack(evalStack, inputs, chromosome)

    def orderStack(self):
        orderStack = self.head + self.tail
        arity = 1
        while arity > 0:
            evalElement = []
            for i in range(arity):
                evalElement.append(orderStack.pop(0))

            arity = 0
            for element in evalElement:
                arity += self.genome.symbols()[element][1]

            self.elementLayers.append(evalElement)

    def evalRecur(self, layerIndex = 0):
        evalStack = []
        for i in range(self.genome.symbols()[self.elementLayers[layerIndex][0]][1]):
            evalStack = evalStack + self.evalRecur(layerIndex + 1)

        evalStack.append(self.elementLayers[layerIndex].pop(0))
        return evalStack

    def evalStack(self, stack, inputs, chromosome):
        returnStack = []
        for symbol in stack:
            if symbol in self.genome.terminals:
                returnStack.append(inputs[symbol])
            elif isNumber(symbol):
                returnStack.append(chromosome.genes[int(symbol)].eval(inputs, chromosome))
            else:
                if self.genome.functions[symbol](returnStack) == 1:
                    return "Error"

        return returnStack.pop()

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False