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
    Ynm_2 = computeYlm(theta1, phi1, lmax=4)
    Ynm_3 = computeYlm(theta2, phi2, lmax=4)
    for nu1 in [0, 2, 4]
        for nu2 in [0, 2, 4]
            for nu3 in [0, 2, 4]
                pre = (
                    (-1)^(nu2 + nu3) /
                    (sqrt(2 * nu3 + 1)) *
                    A(nu2, nu3, nu1, S1, γ1, S2) *
                    A(nu3, S2, γ2, S3)
                )
                for q1 in range(-nu1, nu1, step=2)
                    orientation = B_lpol(nu1, q1, S0, γ0, S1)
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
