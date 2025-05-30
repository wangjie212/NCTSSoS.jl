using DynamicPolynomials
using BenchmarkTools, Profile

@ncpolyvar x y z

mono1 = x^3*y^3
mono2 = x^2*y^2

# no allocation if degrees are not equal
@btime mono1 < mono2

mono3 = x*y^3
# no allocation if variables are the same but time is longer
@btime mono2 < mono3


mono4 = y^3*z
# a ton of allocation if variables are different
# what if we fix the number of variables across all monomials?
@btime mono3 < mono4