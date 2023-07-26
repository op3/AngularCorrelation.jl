module AngularCorrelation

export W, Wcorr, W_coeff, Wcorr_coeff, Parity, EMCharacter, State, Transition, Dipole, Quadrupole, E1, M1, E2, M2

using SphericalHarmonics: Unnormalized, associatedLegendre, sphericalharmonic

include("coefficients.jl")

"""
Applies deorientation and angular distribution coefficients
that result in the final (observed) state after multiple
unobserved transitions.
"""
function final_state_orientation(λ2, cascade...)
    @assert length(cascade) % 2 ≡ 1
    @assert length(cascade) ≥ 3
    res = 1
    for i in 1:2:length(cascade)-4
        Si, γ, Sf = cascade[i:i+2]
        res *= U(λ2, Si, γ, Sf)
    end
    res * A(λ2, cascade[end-2], cascade[end-1], cascade[end])
end

function Wcorr_coeff(
    S0::State, γ0::Transition,
    S1::State, γ1::Transition,
    cascade...)
    coeff = Dict{Vector{Int8},Float64}()
    for λ0 in 0:2:4
        for λ1 in 0:2:4
            for λ2 in 0:2:4
                pre = (
                    (-1)^(λ1 + λ2) /
                    (√(2 * λ2 + 1)) *
                    A(λ1, λ2, λ0, S1, γ1, cascade[begin]) *
                    final_state_orientation(λ2, cascade...)
                )
                for q0 in range(-λ0, λ0, step=2)
                    orientation = B_lpol(q0, λ0, S0, γ0, S1)
                    for q1 in range(-λ1, λ1, step=2)
                        for q2 in range(-λ2, λ2, step=2)
                            val = (
                                pre * orientation *
                                wigner3j(λ2, λ1, λ0, q2, q1, q0)
                            )
                            if val != 0.0
                                if !([λ1, λ2, q1, q2] in keys(coeff))
                                    coeff[[λ1, λ2, q1, q2]] = 0.0
                                end
                                coeff[[λ1, λ2, q1, q2]] += (
                                    pre * orientation *
                                    wigner3j(λ2, λ1, λ0, q2, q1, q0)
                                )
                            end
                        end
                    end
                end
            end
        end
    end
    return coeff
end

function Wcorr(
    theta1::AbstractFloat, phi1::AbstractFloat,
    theta2::AbstractFloat, phi2::AbstractFloat,
    coeff::Dict{Vector{Int8},Float64})
    res = 0.0
    for ((λ1, λ2, q1, q2), c) in coeff
        res += c *
               sphericalharmonic(theta1, phi1, λ1, q1) *
               sphericalharmonic(theta2, phi2, λ2, q2)
    end
    # for λ1 in 0:2:4
    #     for λ2 in 0:2:4
    #         for q1 in range(-λ1, λ1, step=2)
    #             for q2 in range(-λ2, λ2, step=2)
    #                 if coeff[[λ1, λ2, q1, q2]] != 0.0
    #                     res += (
    #                         coeff[[λ1, λ2, q1, q2]] *
    #                         sphericalharmonic(theta1, phi1, λ1, q1) *
    #                         sphericalharmonic(theta2, phi2, λ2, q2)
    #                     )
    #                 end
    #             end
    #         end
    #     end
    # end
    real(res)
end

"""
Angular correlation for a triple γ cascade

The zeroth (exciting) photon of the cascade γ0 is aligned
with the z axis and linearly polarized in the xz plane.
The first (γ1) and last photons are assumed to be observed.
All photons inbetween are unobserved.
This function is normalized to 4π.

Hamilton, eq. 12.204
"""
function Wcorr(
    theta1::AbstractFloat, phi1::AbstractFloat,
    theta2::AbstractFloat, phi2::AbstractFloat,
    S0::State, γ0::Transition,
    S1::State, γ1::Transition,
    cascade...)
    coeff = Wcorr_coeff(S0, γ0, S1, γ1, cascade...)
    Wcorr(theta1, phi1, theta2, phi2, coeff)
end

precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (BigFloat, BigFloat, BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (BigFloat, BigFloat, BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (BigFloat, BigFloat, BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (BigFloat, BigFloat, BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))

function W_coeff(
    S0::State, γ0::Transition,
    cascade...)
    coeff_Pl0 = Dict{Int8,Float64}()
    coeff_Pl2 = Dict{Int8,Float64}()
    for λ in 0:2:4
        fso = final_state_orientation(λ, cascade...)
        coeff_Pl0[λ] = B(λ, S0, γ0, cascade[begin]) * fso
        if λ ≥ 2
            coeff_Pl2[λ] =
                B_lpol(λ, S0, γ0, cascade[begin]) * fso
        end
    end
    return coeff_Pl0, coeff_Pl2
end

"""
Angular distribution for a single emitted γ

The zeroth (exciting) photon of the cascade γ0 is aligned
with the z axis and linearly polarized in the xz plane.
The last photon of the cascade is assumed to be observed.
All photons inbetween are unobserved.
This function is normalized to 4π.

Hamilton, eq. 12.245
"""
function W(
    theta::AbstractFloat, phi::AbstractFloat,
    S0::State, γ0::Transition,
    cascade...)
    res = 0
    coeff_Pl0, coeff_Pl2 = W_coeff(S0, γ0, cascade...)
    c2phi = cos(2 * phi)
    for λ in 0:2:4
        res += coeff_Pl0[λ] * associatedLegendre(theta, l=λ, m=0, norm=Unnormalized())
    end
    for λ in 2:2:4
        res += coeff_Pl2[λ] * associatedLegendre(theta, l=λ, m=2, norm=Unnormalized()) * c2phi
    end
    return res
end

precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State))
precompile(W, (BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (BigFloat, BigFloat, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))

end # module AngularCorrelation
