@enum ComplexKind REAL COMPLEX

struct Variable
    name::String
    kind::ComplexKind

    function Variable(name::AbstractString, kind::ComplexKind=REAL)
        new(convert(String, name), kind)
    end

    Variable(from::Variable, kind::ComplexKind) = new(from.name, kind)
end

function polyarrayvar(complex_kind, prefix, indices...)
    return map(
        i -> Variable("$(prefix)[$(join(i, ","))]", complex_kind),
        Iterators.product(indices...),
    )
end

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

macro ncpolyvar(args...)
    vars, exprs = buildpolyvars(args, REAL)
    # calls the exprs that initializes the variables
    return :($(foldl((x, y) -> :($x; $y), exprs; init=:()));
    $(Expr(:tuple, esc.(vars)...))) # returns the variables
end

# from https://scientificcoder.com/user-defined-show-method-in-julia
# 2-argument show, used by Array show, print(obj) and repr(obj), keep it short
function Base.show(io::IO, obj::Variable)
    print_object(io, obj; multiline=false)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::Variable)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    print_object(io, obj; multiline=multiline)
end

function print_object(io::IO, obj::Variable; multiline::Bool)
    # write something short, or go back to default mode
    multiline ? print(io, "$(obj.name)") : Base.show_default(io, obj)
end
