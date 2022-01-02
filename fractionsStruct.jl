
""" Rational type but using float numbers for numerator and denumerator (specific to Parametric Relaxation) """

# All functions are public

# Frac objects help minimize the number of roundings we do in the parametric linear relaxation.
# Frac is only used in the file parametricLinearRelax which needs exact precision on the computations.

# TYPE

struct Frac
    num::Float64
    den::Float64
    isInteger::Bool
end

# CONSTRUCTOR

function Frac(num,den)
    @assert den != 0 "You can't divide by zero"
    isInteger = false
    if num % 1 == 0.0 && den % 1 == 0.0
        isInteger = true
        gcd_value = gcd(Int(num),Int(den))
        if gcd_value != 1
            num = num/gcd_value
            den = den/gcd_value
        end
    end
    if num < 0 && den < 0
        num = -num
        den = -den
    end
    return Frac(Float64(num),Float64(den),isInteger)
end

# FUNCTIONS AND TESTS

function value_frac(x::Frac)
    return Float64(x.num) / x.den
end

function isZero(x::Frac)
    return x.num == 0
end

function isOne(x::Frac)
    return x.num == x.den
end

function isBetweenZeroAndOne(x::Frac)
    return (x.num > 0 && x.den > 0 && x.num < x.den) || (x.num < 0 && x.den < 0 && x.num > x.den) 
end

# OPERATIONS

function neg_frac(x::Frac)
    return Frac(-x.num,y.den)
end

function inv_frac(x::Frac)
    @assert x.num != 0 "You can't divide by zero"
    return Frac(x.den,x.num)
end

function sub_frac(x::Frac,y::Frac)
    return Frac(x.num * y.den - y.num * x.den,x.den * y.den)
end

function add_frac(x::Frac,y::Frac)
    return Frac(x.num * y.den + y.num * x.den,x.den * y.den)
end

function mult_frac(x::Frac,y::Frac)
    return Frac(x.num*y.num,x.den*y.den)
end

# IMPORT BASE OPERATORS

import Base.:*
import Base.:+
import Base.:-
import Base.:/
import Base.:(==)
import Base.:(!=)
import Base.:(>)
import Base.:(<)
import Base.:(>=)
import Base.:(<=)
import Base.:(isless)

# REWRITE BASE OPERATORS

(-)(x::Frac) = neg_frac(x)
(-)(x::Frac, y::Frac) = sub_frac(x, y)
(+)(x::Frac, y::Frac) = add_frac(x, y)
(*)(x::Frac, y::Frac) = mult_frac(x, y)
(/)(x::Frac, y::Frac) = mult_frac(x, inv_frac(y))

Base.:(==)(x::Frac, y::Frac) = x.num * y.den == x.den * y.num
Base.:(!=)(x::Frac, y::Frac) = x.num * y.den != x.den * y.num
Base.:(>)(x::Frac, y::Frac) = Float64(x.num * y.den - x.den * y.num) / (x.den * y.den) > 0
Base.:(<)(x::Frac, y::Frac) = Float64(x.num * y.den - x.den * y.num) / (x.den * y.den) < 0
Base.:(>=)(x::Frac, y::Frac) = x > y || x == y
Base.:(<=)(x::Frac, y::Frac) = x < y || x == y

function Base.isless(x::Frac, y::Frac)
    return x <= y
end

# SHOW

function Base.show(io::IO, elt::Frac)
    print(io,"$(elt.num) / $(elt.den) ($(value_frac(elt)))")
end

function Base.show(vec::Vector{Frac})
    for elt in vec
        println(elt)
    end
end

function Base.show(vec::Matrix{Frac})
    for elt in vec
        println(elt)
    end
end