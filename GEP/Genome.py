__author__ = 'Coleman'

class Genome():
    def _add(inputs):
        return sum(inputs)
    def _sub(inputs):
        return inputs[0] - inputs[1]
    def _mul(inputs):
        return inputs[0] * inputs[1]
    def _div(inputs):
        if inputs[1] != 0:
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
        self.functions = {"arity": 0}
        self.terminals = []

    def symbols(self):
        return self.functions.update(self.terminals.zip([(0, 0) for i in range(len(self.terminals))]))