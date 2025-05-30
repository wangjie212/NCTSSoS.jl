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

mono5 = y*z

# if the degrees are not the same and variables are different, we don't have memory allocation
@btime mono3 < mono5

# go through the call stack of function

Profile.clear()
@profile for _ in 1:10000000
    mono1 < mono2
end
Profile.print(format=:flat)


mono1 = x^3*y*x
mono2 = x^2*y^2
mono3 = x*y^3
using DynamicPolynomials
@which compare(mono1, mono2)
@which isless(mono1, mono2) 
@which mono1 < mono2
@edit mono1 < mono2

# no allocation if 
@benchmark for _ in 1:1000
    if mono2 < mono3
        @show "hey"
    end
end