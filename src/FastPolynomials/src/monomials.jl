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
        return new(vars, z)
    end
end

function consecutive_unique(vars::Vector{Variable})
    return all(i -> vars[i] != vars[i+1], 1:length(vars)-1)
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
    Monomial(vars, z)

Creates a monomial by converting inputs to proper Vector types.

# Arguments
- `vars`: Collection of variables
- `z`: Collection of exponents

# Returns
- `Monomial`: Monomial with collected variables and exponents
"""
Monomial(vars, z) = Monomial(collect(Variable, vars), collect(Int, z))

Monomial(a::Monomial) = a

Monomial(a::Variable) = Monomial([a], [1])

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
        iszero(z) ? "" : "$(v.name)$(map(c -> exponents[c-'0'+1], string(z)))"
    end
end

function print_object(io::IO, obj::Monomial; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

"""
    Base.hash(m::Monomial, u::UInt)

Computes hash value for a monomial based on its variables and exponents.

# Arguments
- `m::Monomial`: The monomial to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value combining variables and exponents
"""
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
function star(m::Monomial)
    return Monomial(reverse(m.vars), reverse(m.z))
end

"""
    neat_dot(x::Monomial, y::Monomial)

Computes the "neat dot" product of two monomials as star(x) * y.

# Arguments
- `x::Monomial`: First monomial
- `y::Monomial`: Second monomial

# Returns
- `Monomial`: Product of star(x) and y
"""
function neat_dot(x::Monomial, y::Monomial)
    return star(x) * y
end
