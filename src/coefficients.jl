using HalfIntegers: HalfInt

include("wigner.jl")
include("cascade.jl")

"""
ordinary F-coefficient

See Hamilton, eq. 12.168

L: multipolarity
Lp: multipolarity (next-order)

"""
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
@inline function F(nu::Int, nu1::Int, nu0::Int, L::UInt, Lp::UInt, J1::HalfInt, J0::HalfInt)
    (-1)^(Lp + nu1 + nu0 + 1) * sqrt(
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
Angular distribution coefficient

See Hamilton, eq. 12.185.
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
Orientation parameter

See Hamilton, eq. 12.228
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
Orientation parameter

See Hamilton, eq. 12.246
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
        # TODO: Equation in Hamilton?
        if (
            gamma.L + gamma.L >= nu &&
            abs(gamma.L - gamma.L) <= nu &&
            gamma.L + gamma.Lp >= nu &&
            abs(gamma.L - gamma.Lp) < nu &&
            gamma.Lp + gamma.Lp >= nu &&
            abs(gamma.Lp - gamma.Lp) <= nu
        )
            -(1 / 2) *
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
