# Invariants
# ALWAYS make sure z does not contain zeros
# ALWAYS make sure vars does not contain consequtive same variables
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

Monomial(vars, z) = Monomial(collect(Variable, vars), collect(Int, z))
Monomial(a::Monomial) = a
Monomial(a::Variable) = Monomial([a], [1])

degree(m::Monomial) = sum(m.z)
variables(m::Monomial) = unique(m.vars)

# 2-argument show, used by Array show, print(obj) and repr(obj), keep it short
function Base.show(io::IO, obj::Monomial)
    return print_object(io, obj; multiline=false)
end

# the 3-argument show used by display(obj) on the REPL
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

function Base.hash(m::Monomial, u::UInt)
    return hash(m.vars, hash(m.z, u))
end

function Base.one(::Monomial)
    return Monomial(Variable[], Int[])
end

function Base.one(::Type{Monomial})
    return Monomial(Variable[], Int[])
end

function Base.isone(m::Monomial)
    return isempty(m.vars) && isempty(m.z)
end
