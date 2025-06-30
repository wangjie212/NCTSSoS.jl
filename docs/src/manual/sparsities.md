# [Sparsities](@id sparsities)
    The main goal is to reduce the number of monomials used in indexing the moment matrix.

## [Correlative Sparsity](@id correlative-sparsity)

    Correlative sparsity refers to the reduction of monomials that are not independent of each other. For example, in a system with two variables x and y, the monomials x^2 and xy are not independent because they are both proportional to x^2. By eliminating such redundant monomials, we can reduce the size of the moment matrix.

## [Term Sparsity](@id term-sparsity)

    Term sparsity refers to the reduction of monomials that are not present in the system. For example, in a system with two variables x and y, the monomial xy^2 is not present because the system only involves x and y to the first power. By eliminating such missing monomials, we can further reduce the size of the moment matrix.
