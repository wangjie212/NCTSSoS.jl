"""
    Monomial

A monomial represented as a product of variables raised to integer powers.
Automatically consolidates consecutive identical variables and removes zero exponents.

# Fields
- `vars::Vector{Variable}`: Variables in the monomial
- `z::Vector{Int}`: Corresponding non-zero exponents

# Invariants
- `z` does not contain zeros
- `vars` does not contain consecutive same variables
- Length of `vars` equals length of `z`
"""
struct Monomial
    vars::Vector{Variable}
    z::Vector{Int}
    function Monomial(vars::Vector{Variable}, z::Vector{Int})
        @assert length(vars) == length(z) "There should be as many variables as exponents, got $(length(vars)) and $(length(z))"
        @assert consecutive_unique(vars) "Variables should be consecutive unique, got $(vars)"
        @assert all(x -> !iszero(x), z) "Exponents should not be zero"
        return new(vars, z)
    end
end

function consecutive_unique(vars::Vector{Variable})
    return all(i -> vars[i] != vars[i + 1], 1:(length(vars) - 1))
end

function monomial(vars::Vector{Variable}, z::Vector{Int})
    if length(vars) != length(z)
        throw(ArgumentError("There should be as many variables as exponents"))
    end
    nz_idcs = filter(x -> !iszero(z[x]), 1:length(z))
    isempty(nz_idcs) && return Monomial(Variable[], Int[])

    nonconseq_rep_vars = [vars[nz_idcs[1]]]
    nonconseq_rep_z = [z[nz_idcs[1]]]
    for nz_idx in view(nz_idcs, 2:length(nz_idcs))
        if nonconseq_rep_vars[end] == vars[nz_idx]
            nonconseq_rep_z[end] += z[nz_idx]
            continue
        else
            push!(nonconseq_rep_vars, vars[nz_idx])
            push!(nonconseq_rep_z, z[nz_idx])
        end
    end

    return Monomial(nonconseq_rep_vars, nonconseq_rep_z)
end

"""
    monomial(vars, z)

Creates a monomial by converting inputs to proper Vector types.

# Arguments
- `vars`: Collection of variables
- `z`: Collection of exponents

# Returns
- `Monomial`: Monomial with collected variables and exponents
"""
monomial(vars, z) = monomial(collect(Variable, vars), collect(Int, z))

monomial(a::Monomial) = a

monomial(a::Variable) = Monomial([a], [1])

degree(m::Monomial) = sum(m.z)

"""
    variables(m::Monomial)

Extracts all unique variables appearing in a monomial.

# Arguments
- `m::Monomial`: The monomial

# Returns
- `Vector{Variable}`: Unique variables in the monomial
"""
variables(m::Monomial) = unique(m.vars)

function Base.show(io::IO, obj::Monomial)
    return print_object(io, obj; multiline=true)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::Monomial)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

function Base.string(obj::Monomial)
    # https://stackoverflow.com/a/70451947
    exponents = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
    isempty(obj.vars) && return "1"
    return mapreduce(*, zip(obj.vars, obj.z); init="") do (v, z)
        iszero(z) ? "" : "$(string(v))$(map(c -> exponents[c-'0'+1], string(z)))"
    end
end

function print_object(io::IO, obj::Monomial; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

function Base.hash(m::Monomial, u::UInt)
    return hash(m.vars, hash(m.z, u))
end

Base.one(_::Monomial) = Monomial(Variable[], Int[])
Base.one(::Type{Monomial}) = Monomial(Variable[], Int[])
Base.zero(a::Monomial) = Polynomial([0.0], [one(a)])

function Base.isone(m::Monomial)
    return isempty(m.vars) && isempty(m.z)
end

Base.convert(::Type{Monomial}, a::Variable) = Monomial([a], [1])

"""
    star(m::Monomial)

Computes the adjoint (star) of a monomial by reversing variable order and exponents.

# Arguments
- `m::Monomial`: The monomial to compute the adjoint of

# Returns
- `Monomial`: Adjoint monomial with reversed variables and exponents
"""
star(m::Monomial) = Monomial(reverse(m.vars), reverse(m.z))

"""
    _comm(mono::Monomial, comm_gps::Vector{Vector{Variable}})

Projects a monomial onto commutative groups of variables while maintaining the
order of variables within each group.

We kept it as separate groups because symmetric canonicalize over the product of
groups vs product of symmetric canonicalize of each group is different


# Arguments
- `mono::Monomial`: The monomial to project
- `comm_gps::Vector{Vector{Variable}}`: Vector of sets defining commutative variable groups

# Returns
- `Vector{Monomial}`: Projections of the monomial onto each commutative group

# Example
```jldoctest; setup=:(using NCTSSoS.FastPolynomials; using NCTSSoS.FastPolynomials: _comm)
julia> @ncpolyvar x y; comm_gps = [[x], [y]]
2-element Vector{Vector{Variable}}:
 [x]
 [y]

julia> mono1 = x*y*x*y
x¹y¹x¹y¹

julia> _comm(mono1, comm_gps)
2-element Vector{Monomial}:
 x²
 y²
```
"""
function _comm(mono::Monomial, comm_gps::Vector{Vector{Variable}})
    map(comm_gps) do vars
        result = one(Monomial)
        for (var, expon) in zip(mono.vars, mono.z)
            if var in vars
                _mul_var!(result, var, expon, false)
            end
        end
        result
    end
end

# multiply a variable to a monomial
@inline function _mul_var!(result::Monomial, var::Variable, expo::Int, is_unipotent::Bool)
    # Q: do we need to consider commutative case?
    if is_unipotent
        iseven(expo) && return result
        if length(result.vars) == 0 || var != result.vars[end]   # new variable
            push!(result.vars, var)
            push!(result.z, 1)
        else
            pop!(result.vars)
            pop!(result.z)
        end
    else
        iszero(expo) && return result
        if length(result.vars) == 0 || var != result.vars[end]   # new variable
            push!(result.vars, var)
            push!(result.z, expo)
        else
            result.z[end] += expo
        end
    end
    return result
end

"""
    _unipotent(mono::Monomial)

Applies unipotent transformation to a monomial by reducing exponents modulo 2 iteratively.

# Arguments
- `mono::Monomial`: The monomial to transform

# Returns
- `Monomial`: Unipotent form of the monomial with all exponents reduced to 0 or 1
"""
function _unipotent(mono::Monomial)
    isempty(mono.vars) && return mono
    result = one(Monomial)
    for (var, expo) in zip(mono.vars, mono.z)
        _mul_var!(result, var, expo, true)
    end
    return result
end

"""
    _projective(mono::Monomial)

Applies projective transformation to a monomial by setting all non-zero exponents to 1.

# Arguments
- `mono::Monomial`: The monomial to transform

# Returns
- `Monomial`: Projective form with all non-zero exponents set to 1
"""
function _projective(mono::Monomial)
    prod(zip(mono.vars, mono.z); init=one(Monomial)) do (var, expo)
        var^(iszero(expo) ? expo : one(expo))
    end
end

expval(m::Monomial) = m
