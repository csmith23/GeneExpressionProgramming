# Gene Expression Programming
A Gene Expression Programming Implementation in Python

Mk II v1 is the only completed version

To use, give a list of lengths of the expression genes sequences,
a list of lengths of the homeotic genes sequences,
and a list of one or two genomes, which will be used either for both
groups of genes sequences, or the first for the set of expression sequences
and the second for the homeotic sequences.

To construct a genome, instatiate the class and then add a series of
characters, functions, and arity, or a single terminal. The function should accept a list of
inputs that are all numerical values, and use them to compute.
The ordering of the inputs is based on the order which the node recieved the value,
usually from its leftmost child node to its rightmost child node, if the tree were to be drawn out.

To run, give a fitness function, which will be used to compute the fitness of each organism,
in the form of a chromosome, each generation.
It should have all of the inputs and desired outputs allready builtin, or if being used in an
unsupervised setting, be able to access the members of the organism that it requires to properly rate it.
It will be called every round, being passed an organism, and if the organism requires other members to
be tested properly, it should be of a subclass of the Chromosome class.
