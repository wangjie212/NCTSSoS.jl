struct Monomial
    vars::Vector{Variable}
    z::Vector{Int}

    function Monomial(vars::Vector{Variable}, z::Vector{Int})
        if length(vars) != length(z)
            throw(ArgumentError("There should be as many variables as exponents"))
        end
        return new(vars, z)
    end
end

Monomial(vars, z) = Monomial(collect(Variable, vars), collect(Int, z))

# 2-argument show, used by Array show, print(obj) and repr(obj), keep it short
function Base.show(io::IO, obj::Monomial)
    print_object(io, obj; multiline=false)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::Monomial)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    print_object(io, obj; multiline=multiline)
end

function print_object(io::IO, obj::Monomial; multiline::Bool)
    # https://stackoverflow.com/a/70451947
    exponents = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
    return if multiline
        print(
            io,
            mapreduce(*, zip(obj.vars, obj.z); init="") do (v, z)
                iszero(z) ? "" : "$(v.name)$(map(c -> exponents[c-'0'+1], string(z)))"
            end,
        )
    else
        Base.show_default(io, obj)
    end
end
