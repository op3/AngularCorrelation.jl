module AngularCorrelation

export W, Wcorr, W_sample, W_sample_up_to, Wcorr_sample, Wcorr_sample_mt, Wcorr_sample_up_to, W_coeff, Wcorr_coeff
export Parity, EMCharacter, State, Transition, Dipole, Quadrupole, E1, M1, E2, M2

using SphericalHarmonics: Unnormalized, associatedLegendre, sphericalharmonic
using Base.Threads

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

struct CoeffCorr{T<:Real}
    c::Dict{NTuple{4,Int8},T}
end

function Wcorr_coeff(
    S0::State, γ0::Transition,
    S1::State, γ1::Transition,
    cascade...)
    check_cascade(S0, γ0, S1, γ1, cascade...)
    coeff = CoeffCorr{Float64}(Dict{NTuple{4,Int8},Float64}())
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
                            if val ≠ 0.0
                                if !((λ1, λ2, q1, q2) in keys(coeff.c))
                                    coeff.c[(λ1, λ2, q1, q2)] = 0.0
                                end
                                coeff.c[(λ1, λ2, q1, q2)] += (
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

@inline function Wcorr(
    theta1::T, phi1::T,
    theta2::T, phi2::T,
    coeff::CoeffCorr{Float64}) where {T<:Real}
    res = 0.0
    for ((λ1, λ2, q1, q2), c) in coeff.c
        res += c *
               sphericalharmonic(theta1, phi1, λ1, q1) *
               sphericalharmonic(theta2, phi2, λ2, q2)
    end
    real(res)
end

function Wcorr(
    theta1::Vector{T}, phi1::Vector{T},
    theta2::Vector{T}, phi2::Vector{T},
    coeff::CoeffCorr{Float64}) where {T<:Real}
    Wcorr.(theta1, phi1, theta2, phi2, Ref(coeff))
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
@inline function Wcorr(
    theta1::T, phi1::T,
    theta2::T, phi2::T,
    S0::State, γ0::Transition,
    S1::State, γ1::Transition,
    cascade...) where {T<:Real}
    coeff = Wcorr_coeff(S0, γ0, S1, γ1, cascade...)
    Wcorr(theta1, phi1, theta2, phi2, coeff)
end

function Wcorr(
    theta1::Vector{T}, phi1::Vector{T},
    theta2::Vector{T}, phi2::Vector{T},
    cascade...) where {T<:Real}
    Wcorr.(theta1, phi1, theta2, phi2, cascade...)
end

precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Float64, Float64, Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(Wcorr, (Float64, Float64, Float64, Float64, CoeffCorr{Float64}))
precompile(Wcorr, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, CoeffCorr{Float64}))

function W_coeff(
    S0::State, γ0::Transition,
    cascade...)
    check_cascade(S0, γ0, cascade...)
    coeff = zeros(4)
    for λ in [2, 4]
        fso = final_state_orientation(λ, cascade...)
        coeff[λ-1] = B(λ, S0, γ0, cascade[begin]) * fso
        coeff[λ] = B_lpol(λ, S0, γ0, cascade[begin]) * fso
    end
    return coeff
end

function W(
    theta::T, phi::T,
    coeff::Vector{T}) where {T<:Real}
    c2phi = cos(2 * phi)
    (
        1.0 +
        coeff[1] * associatedLegendre(theta, l=2, m=0, norm=Unnormalized()) +
        coeff[2] * associatedLegendre(theta, l=2, m=2, norm=Unnormalized()) * c2phi +
        coeff[3] * associatedLegendre(theta, l=4, m=0, norm=Unnormalized()) +
        coeff[4] * associatedLegendre(theta, l=4, m=2, norm=Unnormalized()) * c2phi
    )
end

function W(
    theta::Vector{T}, phi::Vector{T},
    coeff::Vector{T}) where {T<:Real}
    W.(theta, phi, Ref(coeff))
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
    theta::T, phi::T,
    S0::State, γ0::Transition,
    cascade...) where {T<:Real}
    coeff = W_coeff(S0, γ0, cascade...)
    W(theta, phi, coeff)
end

precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Float64, Float64, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Vector{Float64}, Vector{Float64}, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State, Transition, State))
precompile(W, (Float64, Float64, Float64, Float64, Vector{Float64}))
precompile(W, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}))

function W_estimate_max(coeff)
    theta = 0:π/6:π
    phi = 0:π/6:2π
    vals = W.(theta, phi', Ref(coeff))
    maximum(vals)
end

function Wcorr_estimate_max(coeff)
    theta1 = reshape(Vector(0:π/6:π), (7, 1, 1, 1))
    phi1 = reshape(Vector(0:π/6:2π), (1, 13, 1, 1))
    theta2 = reshape(Vector(0:π/6:π), (1, 1, 7, 1))
    phi2 = reshape(Vector(0:π/6:2π), (1, 1, 1, 13))
    vals = Wcorr.(theta1, phi1, theta2, phi2, Ref(coeff))
    maximum(vals)
end

include("rand_sphere.jl")

function W_sample_up_to(n::Int, coeff::Vector{T}) where {T<:Real}
    Wmax = W_estimate_max(coeff)
    theta, phi = rand_sphere(n)
    values = W.(theta, phi, Ref(coeff))
    sample = rand(n)
    @assert all(values .< Wmax)
    selection = sample .< values
    theta[selection], phi[selection]
end

function W_sample(n::Int, coeff::Vector{T}) where {T<:Real}
    theta = Vector{T}(undef, n)
    phi = Vector{T}(undef, n)
    remaining = n
    sampled = 0

    while remaining > 0
        collected = n - remaining
        efficiency = if sampled ≠ 0
            collected / sampled
        else
            1.0
        end
        sample_next = max(Int64(round(remaining / efficiency * 1.02)), 64)
        nt, np = W_sample_up_to(sample_next, coeff)
        sampled += sample_next
        size = length(nt)

        use = min(remaining, size)
        theta[begin+n-remaining:n-remaining+use] = nt[begin:use]
        phi[begin+n-remaining:n-remaining+use] = np[begin:use]
        remaining -= size
    end

    theta, phi
end

function Wcorr_sample_up_to(n::Int, coeff::CoeffCorr{T}) where {T<:Real}
    Wmax = Wcorr_estimate_max(coeff)
    theta1, phi1 = rand_sphere(n)
    theta2, phi2 = rand_sphere(n)
    values = Wcorr.(theta1, phi1, theta2, phi2, Ref(coeff))
    sample = rand(n)
    @assert all(values .< Wmax)
    selection = sample .< values
    theta1[selection], phi1[selection], theta2[selection], phi2[selection]
end

function Wcorr_sample(n::Int, coeff::CoeffCorr{T}) where {T<:Real}
    theta1 = Vector{T}(undef, n)
    phi1 = Vector{T}(undef, n)
    theta2 = Vector{T}(undef, n)
    phi2 = Vector{T}(undef, n)
    remaining = n
    sampled = 0

    while remaining > 0
        collected = n - remaining
        efficiency = if sampled ≠ 0
            collected / sampled
        else
            1.0
        end
        sample_next = max(Int64(round(remaining / efficiency * 1.02)), 64)
        nt1, np1, nt2, np2 = Wcorr_sample_up_to(sample_next, coeff)
        sampled += sample_next
        size = length(nt1)

        use = min(remaining, size)
        theta1[begin+n-remaining:n-remaining+use] = nt1[begin:use]
        phi1[begin+n-remaining:n-remaining+use] = np1[begin:use]
        theta2[begin+n-remaining:n-remaining+use] = nt2[begin:use]
        phi2[begin+n-remaining:n-remaining+use] = np2[begin:use]
        remaining -= size
    end

    theta1, phi1, theta2, phi2
end

function Wcorr_sample_mt(n::Int, threads::Int, coeff::CoeffCorr{T}) where {T<:Real}
    n_per_thread = Int(n / threads)
    theta1 = Vector{T}(undef, n)
    phi1 = Vector{T}(undef, n)
    theta2 = Vector{T}(undef, n)
    phi2 = Vector{T}(undef, n)

    @threads for i = 1:threads
        tt1, tp1, tt2, tp2 = Wcorr_sample(n_per_thread, coeff)
        theta1[begin+(i-1)*n_per_thread:i*n_per_thread] = tt1
        phi1[begin+(i-1)*n_per_thread:i*n_per_thread] = tp1
        theta2[begin+(i-1)*n_per_thread:i*n_per_thread] = tt2
        phi2[begin+(i-1)*n_per_thread:i*n_per_thread] = tp2
    end
    theta1, phi1, theta2, phi2
end

precompile(W_sample_up_to, (Int64, Vector{Float64}))
precompile(Wcorr_sample_up_to, (Int64, CoeffCorr{Float64}))
precompile(W_sample, (Int64, Vector{Float64}))
precompile(Wcorr_sample, (Int64, CoeffCorr{Float64}))

end # module AngularCorrelation
