# Invariants
# ALWAYS
struct Monomial
    vars::Vector{Variable}
    z::Vector{Int}

    function Monomial(vars::Vector{Variable}, z::Vector{Int})
        if length(vars) != length(z)
            throw(ArgumentError("There should be as many variables as exponents"))
        end
        nz_idcs = filter(x -> !iszero(z[x]), 1:length(z))
        return new(vars[nz_idcs], z[nz_idcs])
    end
end

Monomial(vars, z) = Monomial(collect(Variable, vars), collect(Int, z))

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
