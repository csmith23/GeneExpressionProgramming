__author__ = 'Coleman'

import random

class Gene():
    def __init__(self, genome, homeotic):
        self.genome = genome
        self.head = []
        self.tail = []
        self.homeotic = homeotic

    def initRand(self, length, numGenes):
        arity = self.genome.functions["arity"]
        functions = list(self.genome.functions.keys())
        functions.pop(functions.index("arity"))
        if self.homeotic:
            terminals = list(range(numGenes))
            if numGenes > list(self.genome.genicTerminals)[-1]:
                self.genome.genicTerminals = range(numGenes)
        else:
            terminals = self.genome.terminals

        self.head = [random.choice(functions + terminals) for i in range(length)]
        self.tail = [random.choice(terminals) for i in range(length * (arity - 1) + 1)]
        self.head[0] = random.choice(functions)
        return self

    def replicate(self):
        newGene = Gene(self.genome, self.homeotic)
        newGene.head = self.head[:]
        newGene.tail = self.tail[:]
        return newGene

    def eval(self, inputs):
        elementLayers = self.orderStack()
        evalStack = self.evalRecur(elementLayers)
        return self.evalStack(evalStack, inputs)

    def orderStack(self):
        elementLayers = []
        orderStack = self.head + self.tail
        arity = 1
        while arity > 0:
            evalElement = []
            for i in range(arity):
                evalElement.append(orderStack.pop(0))

            arity = 0
            for element in evalElement:
                arity += self.genome.symbols()[element]

            elementLayers.append(evalElement)

        return elementLayers

    def evalRecur(self, elementLayers, layerIndex = 0):
        evalStack = []
        for i in range(self.genome.symbols()[elementLayers[layerIndex][0]]):
            evalStack = evalStack + self.evalRecur(elementLayers, layerIndex + 1)

        evalStack.append(elementLayers[layerIndex].pop(0))
        return evalStack

    def evalStack(self, stack, inputs):
        returnStack = []
        for symbol in stack:
            if symbol in inputs.keys():
                returnStack.append(inputs[symbol])
            else:
                error = self.genome.functions[symbol][0](returnStack)
                if error == 1:
                    return "Error"

        return returnStack[0]