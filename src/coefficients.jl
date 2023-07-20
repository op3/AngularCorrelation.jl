using HalfIntegers: HalfInt

include("wigner.jl")
include("cascade.jl")

"""
ordinary F-coefficient

See Hamilton, eq. 12.168

L: multipolarity
Lp: multipolarity (next-order)

"""
@inline function F(λ::Int, L::Int, Lp::Int, J2::HalfInt, J1::HalfInt)
    if (λ ≡ 0)
        # Simple case: Hamilton, eq. 12.169
        (L ≡ Lp)
    elseif (
        L + Lp ≥ λ &&
        abs(L - Lp) ≤ λ &&
        Lp + J1 ≥ J2 &&
        abs(Lp - J1) ≤ J2 &&
        L + J2 ≥ J1 &&
        abs(L - J2) ≤ J1 &&
        abs(λ - J1) ≤ J1
    )
        (
            (-1)^(J1 + J2 - 1) *
            sqrt((2 * λ + 1) * (2 * L + 1) * (2 * Lp + 1) * (2 * J1 + 1)) *
            wigner3j(L, Lp, λ, 1, -1, 0) *
            wigner6j(L, Lp, λ, J1, J1, J2)
        )
    else
        0
    end
end

"""
generalized F-coefficient

Hamilton, eq 12.163
"""
@inline function F(λ::Int, λ1::Int, λ0::Int, L::Int, Lp::Int, J1::HalfInt, J0::HalfInt)
    if (λ ≡ 0 && λ1 ≡ 0 && λ0 ≡ 0)
        (L ≡ Lp)
    elseif (L + Lp ≥ λ && abs(L - Lp) ≤ λ)
        (-1)^(Lp + λ1 + λ0 + 1) * sqrt(
            (2 * J0 + 1) *
            (2 * J1 + 1) *
            (2 * L + 1) *
            (2 * Lp + 1) *
            (2 * λ0 + 1) *
            (2 * λ1 + 1) *
            (2 * λ + 1)
        ) * wigner3j(L, Lp, λ, 1, -1, 0) *
        wigner9j(J1, L, J0, J1, Lp, J0, λ1, λ, λ0)
    else
        0
    end
end

"""
Angular distribution coefficient

See Hamilton, eq. 12.185.
"""
@inline function A(λ::Int, state_from::State, γ::Transition, state_to::State)
    (
        F(λ, γ.L, γ.L, state_to.J, state_from.J) +
        F(λ, γ.L, γ.Lp, state_to.J, state_from.J) * 2 * γ.delta +
        F(λ, γ.Lp, γ.Lp, state_to.J, state_from.J) * γ.delta^2
    ) / (1 + γ.delta^2)
end

"""
Generalized directional distribution coefficient

See Hamilton, eq. 12.205

Only for mixture of two multipolarities
"""
@inline function A(λ1::Int, λ2::Int, λ0::Int, state_from::State, γ::Transition, state_to::State)
    (
        F(λ1, λ2, λ0, γ.L, γ.L, state_to.J, state_from.J) +
        F(λ1, λ2, λ0, γ.L, γ.Lp, state_to.J, state_from.J) *
        2 * γ.delta +
        F(λ1, λ2, λ0, γ.Lp, γ.Lp, state_to.J, state_from.J) *
        γ.delta^2
    ) / (1 + γ.delta^2)
end

"""
Orientation parameter

See Hamilton, eq. 12.228
For unpolarized radiation!
"""
@inline function B(λ::Int, state_from::State, γ::Transition, state_to::State)
    (
        F(λ, γ.L, γ.L, state_from.J, state_to.J) +
        F(λ, γ.L, γ.Lp, state_from.J, state_to.J) *
        (-1)^(γ.L + γ.Lp) * 2 * γ.delta +
        F(λ, γ.Lp, γ.Lp, state_from.J, state_to.J) *
        γ.delta^2
    ) / (1 + γ.delta^2)
end

"""
Orientation parameter

See Hamilton, eq. 12.246
For linearly-polarized radiation!
"""
@inline function B_lpol(λ::Int, state_from::State, γ::Transition, state_to::State)
    Int(γ.em_charp) * (
        F(λ, γ.L, γ.L, state_from.J, state_to.J) *
        kappa(λ, γ.L, γ.L) +
        F(λ, γ.L, γ.Lp, state_from.J, state_to.J) *
        kappa(λ, γ.L, γ.Lp) *
        (-1)^(γ.L + γ.Lp) * 2 * γ.delta +
        F(λ, γ.Lp, γ.Lp, state_from.J, state_to.J) *
        kappa(λ, γ.Lp, γ.Lp) *
        γ.delta^2
    ) / (1 + γ.delta^2)
end

"""
Hamilton, eqs. 12.232 to 12.234
For linearly-polarized radiation!
"""
@inline function B_lpol(q::Int, λ::Int, state_from::State, γ::Transition, state_to::State)
    if q ≡ 0
        B(λ, state_from, γ, state_to)
    elseif abs(q) ≡ 2
        if (
            γ.L + γ.L ≥ λ &&
            abs(γ.L - γ.L) ≤ λ &&
            γ.L + γ.Lp ≥ λ &&
            abs(γ.L - γ.Lp) < λ &&
            γ.Lp + γ.Lp ≥ λ &&
            abs(γ.Lp - γ.Lp) ≤ λ
        )
            -(1 / 2) *
            (
                Int(γ.em_char) *
                F(λ, γ.L, γ.L, state_from.J, state_to.J) *
                wigner3j(γ.L, γ.L, λ, 1, 1, -2) /
                wigner3j(γ.L, γ.L, λ, 1, -1, 0)
                -
                Int(γ.em_charp) *
                F(λ, γ.L, γ.Lp, state_from.J, state_to.J) *
                (-1)^(γ.L + γ.Lp) * 2 * γ.delta *
                wigner3j(γ.L, γ.Lp, λ, 1, 1, -2) /
                wigner3j(γ.L, γ.Lp, λ, 1, -1, 0)
                +
                Int(γ.em_charp) *
                F(λ, γ.Lp, γ.Lp, state_from.J, state_to.J) *
                γ.delta^2 *
                wigner3j(γ.Lp, γ.Lp, λ, 1, 1, -2) /
                wigner3j(γ.Lp, γ.Lp, λ, 1, -1, 0)
            ) / (1 + γ.delta^2)
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
@inline function kappa(λ::Int, L::Int, Lp::Int)
    -sqrt(factorial(λ - 2) / factorial(λ + 2)) *
    wigner3j(L, Lp, λ, 1, 1, -2) /
    wigner3j(L, Lp, λ, 1, -1, 0)
end
