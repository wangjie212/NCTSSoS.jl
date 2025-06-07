# WARNING:
# ALWAYS needs to guarantee `monos` are sorted
# ALWAYS needs to guarantee `coeffs` are non-zero
struct Polynomial{T}
    coeffs::Vector{T}
    # perhaps make it into an ordered set?
    monos::Vector{Monomial}

    function Polynomial(a::Vector{T}, x::Vector{Monomial}) where {T}
        length(a) == length(x) ||
            throw(ArgumentError("There should be as many coefficient than monomials"))
        nz_idx = findall(!iszero, a)
        sort!(nz_idx; by=idx -> x[idx])
        return new{T}(a[nz_idx], x[nz_idx])
    end
end

variables(p::Polynomial) = union(variables.(p.monos)...)

function Base.show(io::IO, obj::Polynomial)
    return print_object(io, obj; multiline=false)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::Polynomial)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

function print_object(io::IO, obj::Polynomial; multiline::Bool)
    if multiline
        return print(io, join(
            map(zip(obj.coeffs, obj.monos)) do (coef, mono)
                "$(coef) * $(string(mono))"
            end,
            " + ",
        ))
    else
        return Base.show_default(io, obj)
    end
end

function Base.hash(p::Polynomial, u::UInt)
    return hash(p.coeffs, hash(p.monos, u))
end
