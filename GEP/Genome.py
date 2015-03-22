__author__ = 'Coleman'

class Genome():
    def _add(inputs):
        inputs.append(inputs.pop(-2) + inputs.pop(-1))
    def _sub(inputs):
        inputs.append(inputs.pop(-2) - inputs.pop(-1))
    def _mul(inputs):
        inputs.append(inputs.pop(-2) * inputs.pop(-1))
    def _div(inputs):
        if inputs[1] != 0:
            inputs.append(inputs.pop(-2) / inputs.pop(-1))
        else:
            return 1

    def _NOT(inputs):
        inputs.append(not inputs.pop(-1))
    def _AND(inputs):
        inputs.append(inputs.pop(-2) and inputs.pop(-1))
    def _OR(inputs):
        inputs.append(inputs.pop(-2) or inputs.pop(-1))
    def _XOR(inputs):
        inputs.append(inputs.pop(-2) is not inputs.pop(-1))
    def _XNOR(inputs):
        inputs.append(inputs.pop(-2) is inputs.pop(-1))

    def _exp(inputs):
        inputs.append(inputs.pop(-1) ** 2)
    def _sqrt(inputs):
        inputs.append(inputs.pop(-1) ** 0.5)

    ARITHMETIC_SET = {"+": (_add, 2), "-": (_sub, 2), "*": (_mul, 2), "/": (_div, 2), "arity": 2}
    BOOLEAN_SET = {"N": (_NOT, 1), "A": (_AND, 2), "O": (_OR, 2), "arity": 2}
    REEDMULLEN_SET = {"A": (_AND, 2), "X": (_XOR, 2), "N": (_NOT, 1), "arity": 2}
    EXPONENTIAL_SET = {"E": (_exp, 1), "Q": (_sqrt, 1), "arity": 1}

    def __init__(self):
        self.functions = {"arity": 0}
        self.terminals = []

    def symbols(self):
        return self.functions.update(self.terminals.zip([(0, 0) for i in range(len(self.terminals))]))