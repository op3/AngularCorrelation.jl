using HalfIntegers: twice
using CGcoefficient: f3j, f6j, f9j

"""
Wrapper for f3j to accept parameters without multiplication by two
"""
@inline function wigner3j(
    j1, j2, j3,
    j4, j5, j6)
    t_j1 = convert(Int, twice(j1))
    t_j2 = convert(Int, twice(j2))
    t_j3 = convert(Int, twice(j3))
    t_j4 = convert(Int, twice(j4))
    t_j5 = convert(Int, twice(j5))
    t_j6 = convert(Int, twice(j6))
    f3j(t_j1, t_j2, t_j3, t_j4, t_j5, t_j6)
end

"""
Wrapper for f6j to accept parameters without multiplication by two
"""
@inline function wigner6j(
    j1, j2, j3,
    j4, j5, j6)
    t_j1 = convert(Int, twice(j1))
    t_j2 = convert(Int, twice(j2))
    t_j3 = convert(Int, twice(j3))
    t_j4 = convert(Int, twice(j4))
    t_j5 = convert(Int, twice(j5))
    t_j6 = convert(Int, twice(j6))
    f6j(t_j1, t_j2, t_j3, t_j4, t_j5, t_j6)
end

"""
Wrapper for f9j to accept parameters without multiplication by two
"""
@inline function wigner9j(
    j1, j2, j3,
    j4, j5, j6,
    j7, j8, j9)
    t_j1 = convert(Int, twice(j1))
    t_j2 = convert(Int, twice(j2))
    t_j3 = convert(Int, twice(j3))
    t_j4 = convert(Int, twice(j4))
    t_j5 = convert(Int, twice(j5))
    t_j6 = convert(Int, twice(j6))
    t_j7 = convert(Int, twice(j7))
    t_j8 = convert(Int, twice(j8))
    t_j9 = convert(Int, twice(j9))
    f9j(t_j1, t_j2, t_j3, t_j4, t_j5, t_j6, t_j7, t_j8, t_j9)
end
