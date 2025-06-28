"""
    ComplexKind

An `Enum` representing whether a variable is complex or real
"""
@enum ComplexKind REAL COMPLEX

"""
    Variable

A `Variable` represents a symbolic variable in a polynomial expression.
It can be either `REAL` or `COMPLEX`.
"""
struct Variable
    name::String
    kind::ComplexKind

    function Variable(name::AbstractString, kind::ComplexKind=REAL)
        return new(convert(String, name), kind)
    end

    Variable(from::Variable, kind::ComplexKind) = new(from.name, kind)
end

"""
    polyarrayvar(complex_kind, prefix, indices...)

Creates an array of variables with indexed names.

# Arguments
- `complex_kind::ComplexKind`: The kind of variables (REAL or COMPLEX)
- `prefix::String`: The base name for the variables
- `indices...`: Variable number of index ranges like `i1, i2, ...`

# Returns
- Array of `Variable` objects with names formatted as `prefix[i1,i2,...]`
"""
function polyarrayvar(complex_kind, prefix, indices...)
    return map(
        i -> Variable("$(prefix)[$(join(i, ","))]", complex_kind),
        Iterators.product(indices...),
    )
end

"""
    buildpolyvar(var, complex_kind)

Builds a polynomial variable declaration from a symbol or expression.

# Arguments
- `var`: Either a Symbol for single variable or Expr for array variables
- `complex_kind::ComplexKind`: The kind of variables (REAL or COMPLEX)

# Returns
- Tuple of (variable_name, expression) for variable creation
"""
function buildpolyvar(var, complex_kind)
    if isa(var, Symbol)
        var, :($(esc(var)) = $Variable($"$var", $complex_kind))
    else
        isa(var, Expr) || error("Expected $var to be a variable name")
        Base.Meta.isexpr(var, :ref) ||
            error("Expected $var to be of the form varname[idxset]")
        (2 ≤ length(var.args)) || error("Expected $var to have at least one index set")
        varname = var.args[1]
        prefix = string(varname)
        varname,
        :(
            $(esc(varname)) = polyarrayvar(
                $(complex_kind), $prefix, $(esc.(var.args[2:end])...)
            )
        )
    end
end

"""
    buildpolyvars(args, complex_kind)

Builds multiple polynomial variable declarations from a collection of arguments.

# Arguments
- `args`: Collection of symbols or expressions for variable names
- `complex_kind::ComplexKind`: The kind of variables (REAL or COMPLEX)

# Returns
- Tuple of (variable_names_array, expressions_array) for variable creation
"""
function buildpolyvars(args, complex_kind)
    vars = Symbol[]
    exprs = []
    for arg in args
        var, expr = buildpolyvar(arg, complex_kind)
        push!(vars, var)
        push!(exprs, expr)
    end
    return vars, exprs
end

"""
    @ncpolyvar(args...)

Macro to create non-commutative polynomial variables of REAL kind.

# Arguments
- `args...`: Variable arguments specifying variable names (symbols or indexed expressions)

# Returns
- Tuple of created Variable objects

# Example
```jldoctest; setup=:(using NCTSSoS.FastPolynomials)
julia> @ncpolyvar x y z     # Creates three real variables
(x, y, z)

julia> @ncpolyvar u[1:3]     # Creates array of variables u[1], u[2], u[3]
(Variable[u[1], u[2], u[3]],)
```

"""
macro ncpolyvar(args...)
    vars, exprs = buildpolyvars(args, REAL)
    # calls the exprs that initializes the variables
    return :($(foldl((x, y) -> :($x; $y), exprs; init=:()));
    $(Expr(:tuple, esc.(vars)...))) # returns the variables
end

# from https://scientificcoder.com/user-defined-show-method-in-julia
# 2-argument show, used by Array show, print(obj) and repr(obj), keep it short
function Base.show(io::IO, obj::Variable)
    return print_object(io, obj; multiline=true)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::Variable)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

"""
    print_object(io, obj; multiline)

Prints a Variable object to an IO stream.

# Arguments
- `io::IO`: The output stream
- `obj::Variable`: The Variable to print
- `multiline::Bool`: Whether to use multiline format

# Returns
- Nothing (prints to IO stream)
"""
function print_object(io::IO, obj::Variable; multiline::Bool)
    return multiline ? print(io, "$(obj.name)") : Base.show_default(io, obj)
end

"""
    Base.hash(v, u)

Computes hash value for a Variable based on its name.

# Arguments
- `v::Variable`: The Variable to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value of the variable's name
"""
function Base.hash(v::Variable, u::UInt)
    return hash(v.name, u)
end

Base.one(::Variable) = Monomial(Variable[], Int[])

"""
    monomials(vars::Vector{Variable}, cur_d::Int)

Generates all monomials of a specific degree from given variables.

# Arguments
- `vars::Vector{Variable}`: Variables to use in monomials
- `cur_d::Int`: Degree of monomials to generate

# Returns
- `Vector{Monomial}`: All monomials of degree `cur_d` in the given variables
"""
function monomials(vars::Vector{Variable}, ::Val{D}) where {D}
    return map(cartesian_product(vars, Val{D}())) do cur_vars
        monomial(collect(Variable, cur_vars), ones(Int, length(cur_vars)))
    end
end

function cartesian_product(vars::Vector{T}, ::Val{D}) where {T,D}
    iterable = Iterators.product(ntuple(_ -> vars, Val{D}())...)
    a = collect(iterable)
    return reshape(a, :)
end

"""
    get_basis(vars::Vector{Variable}, d::Int)

Generates a sorted basis of all monomials up to a given degree.

# Arguments
- `vars::Vector{Variable}`: Variables to use in the basis
- `d::Int`: Maximum degree of monomials

# Returns
- `Vector{Monomial}`: Sorted basis containing all monomials of degrees `0` through `d`
"""
function get_basis(vars::Vector{Variable}, d::Int)
    vec_of_monos = Vector{Vector{Monomial}}(undef, d + 1)
    for i in 0:d
        vec_of_monos[i + 1] = monomials(vars, Val(i))
    end
    return sort!(reduce(vcat, vec_of_monos))
end


function Base.:(^)(a::Variable, expo::Int)
    @assert expo >= 0 "Exponent must be non-negative."
    return iszero(expo) ? one(a) : Monomial([a], [expo])
end

function Base.string(a::Variable)
    subscripts = ('₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉')
    m = match(r"(?<var>\w+)(?:\[(?<idx>\d+)\])?",a.name)
    isnothing(m[:idx]) ? m[:var] :  m[:var] * map(c->subscripts[c-'0'+1], m[:idx])
end
