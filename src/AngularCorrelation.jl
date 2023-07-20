module AngularCorrelation

export Wcorr, Parity, EMCharacter, State, Transition, Dipole, Quadrupole, E1, M1, E2, M2

using EnumX
using HalfIntegers: half, Half, HalfInt, twice
using RationalRoots: signedroot, RationalRoot
using CGcoefficient: f3j, f6j, f9j
using SphericalHarmonics: computeYlm

@enumx EMCharacter electric = -1 magnetic = 1 unknown = 0

struct Transition
    L::UInt
    em_char::EMCharacter.T
    Lp::UInt
    em_charp::EMCharacter.T
    delta::Real
end

"""
Invert the radiation character

char: original radiation character
"""
function alt_char(char::EMCharacter.T)
    if char == Parity.negative
        Parity.positive
    elseif char == Parity.positive
        Parity.negative
    else
        Parity.unknown
    end
end

E1(delta=0.0) = Transition(HalfInt(1), EMCharacter.electric, HalfInt(2), EMCharacter.magnetic, delta)
M1(delta=0.0) = Transition(HalfInt(1), EMCharacter.magnetic, HalfInt(2), EMCharacter.electric, delta)
E2(delta=0.0) = Transition(HalfInt(2), EMCharacter.electric, HalfInt(3), EMCharacter.magnetic, delta)
M2(delta=0.0) = Transition(HalfInt(2), EMCharacter.magnetic, HalfInt(3), EMCharacter.electric, delta)
Dipole(delta=0.0) = Transition(HalfInt(1), EMCharacter.unknown, HalfInt(2), EMCharacter.unknown, delta)
Quadrupole(delta=0.0) = Transition(HalfInt(2), EMCharacter.unknown, HalfInt(3), EMCharacter.unknown, delta)

Base.broadcastable(x::Transition) = Ref(x)

@enumx Parity negative = -1 positive = 1 unknown = 0

struct State{T}
    J::Half{T}
    parity::Parity.T
end

Base.broadcastable(x::State) = Ref(x)

function State(J::Int, parity::Parity.T)
    State(HalfInt(J), parity)
end

function State(J::Rational, parity::Parity.T)
    State(HalfInt(J), parity)
end

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


"""
ordinary F-coefficient

See Hamilton, eq. 12.168

L: multipolarity
Lp: multipolarity (next-order)

"""
#const lru_F_ord = LRU{Tuple{Int, UInt, UInt, HalfInt, HalfInt}, Float64}(maxsize=20000)
#function F(nu::Int, L::UInt, Lp::UInt, J2::HalfInt, J1::HalfInt)
#    get!(lru_F_ord, (nu, L, Lp, J2, J1)) do
#        _F(nu, L, Lp, J2, J1)
#    end
#end

@inline function F(nu::Int, L::UInt, Lp::UInt, J2::HalfInt, J1::HalfInt)
    if (nu == 0)
        (L == Lp)
    else
        (
            (-1)^(J1 + J2 - 1) *
            sqrt((2 * nu + 1) * (2 * L + 1) * (2 * Lp + 1) * (2 * J1 + 1)) *
            wigner3j(L, Lp, nu, 1, -1, 0) *
            wigner6j(L, Lp, nu, J1, J1, J2)
        )
    end 
end

"""
generalized F-coefficient

Hamilton, eq 12.163
"""
#const lru_F_gen = LRU{Tuple{Int, Int, Int, UInt, UInt, HalfInt, HalfInt}, Float64}(maxsize=20000)
#function F(nu::Int, nu1::Int, nu0::Int, L::UInt, Lp::UInt, J2::HalfInt, J1::HalfInt)
#    get!(lru_F_gen, (nu, nu1, nu0, L, Lp, J2, J1)) do
#        _F(nu, nu1, nu0, L, Lp, J2, J1)
#    end
#end

@inline function F(nu::Int, nu1::Int, nu0::Int, L::UInt, Lp::UInt, J1::HalfInt, J0::HalfInt)
    (-1)^(Lp+nu1+nu0+1) * sqrt(
        (2 * J0 + 1) *
        (2 * J1 + 1) *
        (2 * L + 1) *
        (2 * Lp + 1) *
        (2 * nu0 + 1) *
        (2 * nu1 + 1) *
        (2 * nu + 1)
    ) * wigner3j(L, Lp, nu, 1, -1, 0) *
    wigner9j(J1, L, J0, J1, Lp, J0, nu1, nu, nu0)
end

"""
See Hamilton, eq. 12.185.

Angular distribution coefficient
"""
@inline function A(nu::Int, state_from::State, gamma::Transition, state_to::State)
    (
        F(nu, gamma.L, gamma.L, state_to.J, state_from.J) +
        F(nu, gamma.L, gamma.Lp, state_to.J, state_from.J) * 2 * gamma.delta +
        F(nu, gamma.Lp, gamma.Lp, state_to.J, state_from.J) * gamma.delta^2
    ) / (1 + gamma.delta^2)
end

"""
Generalized directional distribution coefficient

See Hamilton, eq. 12.205

Only for mixture of two multipolarities
"""
@inline function A(nu1::Int, nu2::Int, nu0::Int, state_from::State, gamma::Transition, state_to::State)
    (
        F(nu1, nu2, nu0, gamma.L, gamma.L, state_to.J, state_from.J) +
        F(nu1, nu2, nu0, gamma.L, gamma.Lp, state_to.J, state_from.J) *
        2 * gamma.delta +
        F(nu1, nu2, nu0, gamma.Lp, gamma.Lp, state_to.J, state_from.J) *
        gamma.delta^2
    ) / (1 + gamma.delta^2)
end

"""
See Hamilton, eq. 12.228
Orientation parameter
For unpolarized radiation!
"""
@inline function B(nu::Int, state_from::State, gamma::Transition, state_to::State)
    (
        F(nu, gamma.L, gamma.L, state_from.J, state_to.J) +
        F(nu, gamma.L, gamma.Lp, state_from.J, state_to.J) *
        (-1)^(gamma.L + gamma.Lp) * 2 * gamma.delta +
        F(nu, gamma.Lp, gamma.Lp, state_from.J, state_to.J) *
        gamma.delta^2
    ) / (1 + gamma.delta^2)
end

"""
See Hamilton, eq. 12.246
Orientation parameter
For linearly-polarized radiation!
"""
@inline function B_lpol(nu::Int, state_from::State, gamma::Transition, state_to::State)
    Int(gamma.em_charp) * (
        F(nu, gamma.L, gamma.L, state_from.J, state_to.J) *
        kappa(nu, gamma.L, gamma.L) +
        F(nu, gamma.L, gamma.Lp, state_from.J, state_to.J) *
        kappa(nu, gamma.L, gamma.Lp) *
        (-1)^(gamma.L + gamma.Lp) * 2 * gamma.delta +
        F(nu, gamma.Lp, gamma.Lp, state_from.J, state_to.J) *
        kappa(nu, gamma.Lp, gamma.Lp) *
        gamma.delta^2
    ) / (1 + gamma.delta^2)
end

"""
Hamilton, eqs. 12.232 to 12.234
For linearly-polarized radiation!
"""
@inline function B_lpol(nu::Int, q::Int, state_from::State, gamma::Transition, state_to::State)
    if q == 0
        B(nu, state_from, gamma, state_to)
    elseif abs(q) == 2
        -(1/2) *
        Int(gamma.em_charp) *
        (
            F(nu, gamma.L, gamma.L, state_from.J, state_to.J) *
            wigner3j(gamma.L, gamma.L, nu, 1, 1, -2) /
            wigner3j(gamma.L, gamma.L, nu, 1, -1, 0)
            +
            F(nu, gamma.L, gamma.Lp, state_from.J, state_to.J) *
            (-1)^(gamma.L + gamma.Lp) * 2 * gamma.delta *
            wigner3j(gamma.L, gamma.Lp, nu, 1, 1, -2) /
            wigner3j(gamma.L, gamma.Lp, nu, 1, -1, 0)
            +
            F(nu, gamma.Lp, gamma.Lp, state_from.J, state_to.J) *
            gamma.delta^2 *
            wigner3j(gamma.Lp, gamma.Lp, nu, 1, 1, -2) /
            wigner3j(gamma.Lp, gamma.Lp, nu, 1, -1, 0)
        ) / (1 + gamma.delta^2)
    else
        0
    end
end

"""
See Hamilton, eq. 12.243
"""
@inline function kappa(nu::Int, L::UInt, Lp::UInt)
    -sqrt(factorial(nu - 2) / factorial(nu + 2)) *
    wigner3j(L, Lp, nu, 1, 1, -2) /
    wigner3j(L, Lp, nu, 1, -1, 0)
end

"""
Stub

Hamilton, eq. 12.204

I0 → I1 → I2 → I3
with γ0 → γ1 → γ2

"""
function Wcorr(
        theta1::Float64, phi1::Float64,
        theta2::Float64, phi2::Float64,
        S0::State, g0::Transition,
        S1::State, g1::Transition,
        S2::State, g2::Transition,
        S3::State)
    res = 0
    Ynm_2 = computeYlm(theta1, phi1, lmax=4)
    Ynm_3 = computeYlm(theta2, phi2, lmax=4)
    for nu1 in [0, 2, 4]
        for nu2 in [0, 2, 4]
            for nu3 in [0, 2, 4]
                pre = (
                    (-1)^(nu2 + nu3) /
                    (sqrt(2 * nu3 + 1)) *
                    A(nu2, nu3, nu1, S1, g1, S2) *
                    A(nu3, S2, g2, S3)
                )
                for q1 in range(-nu1, nu1, step=2)
                    orientation = B_lpol(nu1, q1, S0, g0, S1)
                    for q2 in range(-nu2, nu2, step=2)
                        for q3 in range(-nu3, nu3, step=2)
                            res += (
                                pre * orientation *
                                # This is where U-coefficients would go
                                wigner3j(nu3, nu2, nu1, q3, q2, q1) *
                                Ynm_2[(nu2, q2)] *
                                Ynm_3[(nu3, q3)]
                            )
                        end
                    end
                end
            end
        end
    end
    real(res)
end

precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State))

end # module AngularCorrelation
