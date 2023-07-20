module AngularCorrelation

export Wcorr, Parity, EMCharacter, State, Transition, Dipole, Quadrupole, E1, M1, E2, M2

using SphericalHarmonics: computeYlm

include("coefficients.jl")

"""
Angular correlation for a triple γ cascade

S0 → S1 → S2 → S3
with γ0 → γ1 → γ2

The first photon γ0 is assumed to be aligned
with the z axis and linearly polarized in the xz plane

Hamilton, eq. 12.204
"""
function Wcorr(
    theta1::Float64, phi1::Float64,
    theta2::Float64, phi2::Float64,
    S0::State, γ0::Transition,
    S1::State, γ1::Transition,
    S2::State, γ2::Transition,
    S3::State)
    res = 0
    Ynm_1 = computeYlm(theta1, phi1, lmax=4)
    Ynm_2 = computeYlm(theta2, phi2, lmax=4)
    for λ0 in 0:2:4
        for λ1 in 0:2:4
            for λ2 in 0:2:4
                pre = (
                    (-1)^(λ1 + λ2) /
                    (sqrt(2 * λ2 + 1)) *
                    A(λ1, λ2, λ0, S1, γ1, S2) *
                    A(λ2, S2, γ2, S3)
                )
                for q0 in range(-λ0, λ0, step=2)
                    orientation = B_lpol(q0, λ0, S0, γ0, S1)
                    for q1 in range(-λ1, λ1, step=2)
                        for q2 in range(-λ2, λ2, step=2)
                            res += (
                                pre * orientation *
                                # This is where U-coefficients would go
                                wigner3j(λ2, λ1, λ0, q2, q1, q0) *
                                Ynm_1[(λ1, q1)] *
                                Ynm_2[(λ2, q2)]
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
