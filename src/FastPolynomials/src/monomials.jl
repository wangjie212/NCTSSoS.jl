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
        if length(vars) != length(z)
            throw(ArgumentError("There should be as many variables as exponents"))
        end
        nz_idcs = filter(x -> !iszero(z[x]), 1:length(z))
        isempty(nz_idcs) && return new(Variable[], Int[])

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

        return new(nonconseq_rep_vars, nonconseq_rep_z)
    end
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

"""
    Monomial(a::Monomial)

Identity constructor for monomials.

# Arguments
- `a::Monomial`: Monomial to copy

# Returns
- `Monomial`: The same monomial
"""
Monomial(a::Monomial) = a

"""
    Monomial(a::Variable)

Creates a monomial from a single variable with exponent 1.

# Arguments
- `a::Variable`: Variable to convert

# Returns
- `Monomial`: Monomial representing the variable
"""
Monomial(a::Variable) = Monomial([a], [1])

"""
    degree(m::Monomial)

Computes the total degree of a monomial (sum of all exponents).

# Arguments
- `m::Monomial`: The monomial

# Returns
- `Int`: Total degree of the monomial
"""
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

# 2-argument show, used by Array show, print(obj) and repr(obj), keep it short
function Base.show(io::IO, obj::Monomial)
    return print_object(io, obj; multiline=true)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::Monomial)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

"""
    Base.string(obj::Monomial)

Converts a monomial to its string representation using Unicode superscript exponents.

# Arguments
- `obj::Monomial`: The monomial to convert

# Returns
- `String`: String representation with superscript exponents, "1" for empty monomial
"""
function Base.string(obj::Monomial)
    # https://stackoverflow.com/a/70451947
    exponents = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
    isempty(obj.vars) && return "1"
    return mapreduce(*, zip(obj.vars, obj.z); init="") do (v, z)
        iszero(z) ? "" : "$(v.name)$(map(c -> exponents[c-'0'+1], string(z)))"
    end
end

"""
    print_object(io, obj; multiline)

Prints a Monomial object to an IO stream.

# Arguments
- `io::IO`: The output stream
- `obj::Monomial`: The monomial to print
- `multiline::Bool`: Whether to use string representation or default format

# Returns
- Nothing (prints to IO stream)
"""
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

"""
    Base.one(::Monomial)

Returns the multiplicative identity (empty monomial representing 1).

# Arguments
- `::Monomial`: Any monomial instance (type parameter only)

# Returns
- `Monomial`: Empty monomial representing 1
"""
function Base.one(::Monomial)
    return Monomial(Variable[], Int[])
end

"""
    Base.one(::Type{Monomial})

Returns the multiplicative identity for the Monomial type.

# Arguments
- `::Type{Monomial}`: The Monomial type

# Returns
- `Monomial`: Empty monomial representing 1
"""
function Base.one(::Type{Monomial})
    return Monomial(Variable[], Int[])
end

"""
    Base.isone(m::Monomial)

Checks if a monomial represents the multiplicative identity (1).

# Arguments
- `m::Monomial`: The monomial to check

# Returns
- `Bool`: True if monomial is empty (represents 1), false otherwise
"""
function Base.isone(m::Monomial)
    return isempty(m.vars) && isempty(m.z)
end
