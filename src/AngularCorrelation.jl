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

"""
    W(
        theta::T, phi::T,
        coeff::Vector{T}) where {T<:Real}

Calculate angular distribution for a single emitted γ.

The zeroth (exciting) photon of the cascade γ0 is aligned
with the z axis and linearly polarized in the xz plane.
The last photon of the cascade is assumed to be observed.
All photons inbetween are unobserved.
This function is normalized to 4π.

See Hamilton, eq. 12.245

# Arguments
- `theta::T`: azimuthal angle of emitted photon.
- `phi::T`: polar angle of emitted photon.
- `coeff::Vector{T}`: Angular distribution coefficients obtained from `W_coeff()`.
"""
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

"""
    W(
        theta::Vector{T}, phi::Vector{T},
        coeff::Vector{T}) where {T<:Real}

Calculate angular distribution for a single emitted γ.

The zeroth (exciting) photon of the cascade γ0 is aligned
with the z axis and linearly polarized in the xz plane.
The last photon of the cascade is assumed to be observed.
All photons inbetween are unobserved.
This function is normalized to 4π.

See Hamilton, eq. 12.245

# Arguments
- `theta::Vector{T}`: azimuthal angle of emitted photon.
- `phi::Vector{T}`: polar angle of emitted photon.
- `coeff::Vector{T}`: Angular distribution coefficients obtained from `W_coeff()`.
"""
function W(
    theta::Vector{T}, phi::Vector{T},
    coeff::Vector{T}) where {T<:Real}
    W.(theta, phi, Ref(coeff))
end

"""
    W(
        theta::T, phi::T,
        S0::State, γ0::Transition,
        cascade...) where {T<:Real}

Calculate angular distribution for a single emitted γ.

The zeroth (exciting) photon of the cascade γ0 is aligned
with the z axis and linearly polarized in the xz plane.
The last photon of the cascade is assumed to be observed.
All photons inbetween are unobserved.
This function is normalized to 4π.

See Hamilton, eq. 12.245

# Arguments
- `theta::T`: azimuthal angle of emitted photon.
- `phi::T`: polar angle of emitted photon.
- `S0::State`: Initial state before excitation.
- `γ0::Transition`: Photon that performs the excitation
- `cascade...`: Rest of the decay cascade
  (consisting of `Transition`s and `State`s)
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

"""
    W_estimate_max(coeff)

Estimate the maximum of an angular distribution.
This is done by sampling the distribution at angles
that are multiples or π/6.

# Arguments
- `coeff::Vector{T}`: Coefficients for angular distribution obtained
  using `W_coeff()`

# Examples
```julia-repl
julia> coeff = W_coeff(State(0), E1(), State(1), Dipole(), State(0))
julia> AngularCorrelation.W_estimate_max(coeff)
1.5
```
"""
function W_estimate_max(coeff)
    theta = 0:π/12:π/2
    phi = 0:π/12:π/2
    vals = abs.(W.(theta, phi', Ref(coeff)))
    maximum(vals)
end

function Wcorr_estimate_max(coeff)
    thetas = Vector(0:π/12:π/2)
    phis = Vector(0:π/12:π/2)
    theta1 = reshape(thetas, (length(thetas), 1, 1, 1))
    phi1 = reshape(phis, (1, length(phis), 1, 1))
    theta2 = reshape(thetas, (1, 1, length(thetas), 1))
    phi2 = reshape(phis, (1, 1, 1, length(phis)))
    vals = abs.(Wcorr.(theta1, phi1, theta2, phi2, Ref(coeff)))
    # Safty factor of 0.5% because I have no idea how to properly find
    # the real maximum. The maximum deviation observed in practice is
    # 3.4713%, giving some slight margin for further errors.
    maximum(vals) .* 1.005
end

include("rand_sphere.jl")

"""
    W_sample_up_to(n::Integer, coeff::Vector{T}) where {T<:Real}

Obtain sample for angular distribution.
Retults in an array of size(n, 2)
with the first dimension corresponding to the number of non-rejected samples,
and the second dimension corresponding to the coordinates theta and phi.
The number of non-rejected samples that are returned is smaller or equal
to the number of samples `n` drawn initially.
Depending on how strongly correlated the angular distribution is, significantly
less samples can be returned.

# Arguments
- `n::Integer`: Number of samples to draw
- `coeff::Vector{T}`: Coefficients for angular distribution obtained
  using `W_coeff()`

# Examples
```julia-repl
julia> coeff = W_coeff(State(0), E1(), State(1), Dipole(), State(0))
julia> W_sample_up_to(5, coeff)
3×2 Matrix{Float64}:
 2.00991  6.1304
 2.54216  2.32825
 1.02637  0.953504
 ```
"""
function W_sample_up_to(n::Int, coeff::Vector{T}) where {T<:Real}
    Wmax = W_estimate_max(coeff)
    sph = rand_sphere(n)'
    values = W.(sph[:, 1], sph[:, 2], Ref(coeff))
    sample = rand(n)
    @assert all(values .< Wmax)
    selection = sample .< values
    sph[selection, :]
end

"""
    W_sample(n::Integer, coeff::Vector{T}) where {T<:Real}

Obtain sample for angular distribution.
Retults in an array of size(n, 2)
with the first dimension corresponding to the number of samples,
and the second dimension corresponding to the coordinates theta and phi.

# Arguments
- `n::Integer`: Number of samples to draw
- `coeff::Vector{T}`: Coefficients for angular distribution obtained
  using `W_coeff()`

# Examples
```julia-repl
julia> coeff = W_coeff(State(0), E1(), State(1), Dipole(), State(0))
julia> W_sample(5, coeff)
5×2 Matrix{Float64}:
 1.46641   1.04276
 2.71659   4.13127
 2.7962    5.12808
 0.986001  6.00779
 0.735967  0.32747
 ```
"""
function W_sample(n::Integer, coeff::Vector{T}) where {T<:Real}
    res = Array{T}(undef, n, 2)
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
        sample_new = W_sample_up_to(sample_next, coeff)
        sampled += sample_next
        elements = size(sample_new)[1]

        use = min(remaining, elements)
        res[begin+n-remaining:n-remaining+use, :] = sample_new[begin:use, :]
        remaining -= elements
    end
    res
end

"""
    Wcorr_sample_up_to(n::Integer, coeff::CoeffCorr{T}) where {T<:Real}

Obtain sample for angular correlation. 
Results in an array of size (i, 2, 2)
with the first dimension corresponding to the number of non-rejected samples,
the second dimension corresponding to the detected photon,
and the third dimension corresponding to the coordinates theta and phi.
The number of non-rejected samples that are returned is smaller or equal
to the number of samples `n` drawn initially.
Depending on how strongly correlated the angular distribution is, significantly
less samples can be returned.

# Arguments
- `n::Integer`: Number of samples to draw before rejection
- `coeff::CoeffCorr{T}`: Coefficients for angular correlation obtained
  using `Wcorr_coeff()`

# Examples
```julia-repl
julia> coeff = Wcorr_coeff(State(0), E1(), State(1), Dipole(), State(2), Quadrupole(), State(0))
julia> Wcorr_sample_up_to(50, coeff)
3×2×2 Array{Float64, 3}:
[:, :, 1] =
 1.06189   0.435861
 0.541557  1.23968
 2.57618   1.72244

[:, :, 2] =
 2.66601  4.03896
 1.90058  5.88181
 4.72594  4.09244
```
"""
function Wcorr_sample_up_to(n::Integer, coeff::CoeffCorr{T}) where {T<:Real}
    Wmax = Wcorr_estimate_max(coeff)
    # (samples, gamma, coordinate)
    sph = Array{Float64}(undef, n, 2, 2)
    sph[:, 1, :] = rand_sphere(n)'
    sph[:, 2, :] = rand_sphere(n)'
    values = Wcorr.(sph[:, 1, 1], sph[:, 1, 2], sph[:, 2, 1], sph[:, 2, 2], Ref(coeff))
    sample = rand(n)
    @assert all(values .< Wmax)
    selection = sample .< values
    sph[selection, :, :]
end

"""
    Wcorr_sample(n::Integer, coeff::CoeffCorr{T}) where {T<:Real}

Obtain sample for angular correlation. 
Results in an array of size (n, 2, 2)
with the first dimension corresponding to the number of samples,
the second dimension corresponding to the detected photon,
and the third dimension corresponding to the coordinates theta and phi.

# Arguments
- `n::Integer`: Number of samples to draw
- `coeff::CoeffCorr{T}`: Coefficients for angular correlation obtained
  using `Wcorr_coeff()`

# Examples
```julia-repl
julia> coeff = Wcorr_coeff(State(0), E1(), State(1), Dipole(), State(2), Quadrupole(), State(0))
julia> Wcorr_sample(5, coeff)
5×2×2 Array{Float64, 3}:
[:, :, 1] =
 1.22272  2.08489
 1.26879  1.24111
 2.52387  0.457633
 2.18338  1.12996
 2.55434  0.668774

[:, :, 2] =
 4.41387  2.41454
 1.04299  1.37492
 4.25477  3.68128
 1.90813  3.43944
 2.08456  0.260601
```
"""
function Wcorr_sample(n::Integer, coeff::CoeffCorr{T}) where {T<:Real}
    res = Array{Float64}(undef, n, 2, 2)
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
        sample_new = Wcorr_sample_up_to(sample_next, coeff)
        sampled += sample_next
        elements = size(sample_new)[1]

        use = min(remaining, elements)
        res[begin+n-remaining:n-remaining+use, :, :] = sample_new[begin:use, :, :]
        remaining -= elements
    end
    res
end

"""
    Wcorr_sample_mt(n::Integer, coeff::CoeffCorr{T}) where {T<:Real}

Obtain sample for angular correlation. 
Results in an array of size (n, 2, 2)
with the first dimension corresponding to the number of samples,
the second dimension corresponding to the detected photon,
and the third dimension corresponding to the coordinates theta and phi.
The calculation is multi-threaded, set the environment variable
`JULIA_NUM_THREADS` before starting julia or use
the cli argument `-t/--threads`.

# Arguments
- `n::Integer`: Number of samples to draw
- `threads::Integer`: Number of threads to use
- `coeff::CoeffCorr{T}`: Coefficients for angular correlation obtained
  using `Wcorr_coeff()`

# Examples
```julia-repl
julia> coeff = Wcorr_coeff(State(0), E1(), State(1), Dipole(), State(2), Quadrupole(), State(0))
julia> Wcorr_sample_mt(4, 2, coeff)
4×2×2 Array{Float64, 3}:
[:, :, 1] =
 1.7967   1.63842
 1.88768  1.38813
 1.04846  0.863921
 2.16277  1.1874

[:, :, 2] =
 1.20264  2.84755
 3.24036  3.75469
 4.66713  3.03743
 3.75127  2.38884
```
"""
function Wcorr_sample_mt(n::Integer, threads::Integer, coeff::CoeffCorr{T}) where {T<:Real}
    n_per_thread = Int(n / threads)
    res = Array{Float64}(undef, n, 2, 2)

    @threads for i = 1:threads
        res[begin+(i-1)*n_per_thread:i*n_per_thread, :, :] = Wcorr_sample(n_per_thread, coeff)
    end
    res
end

precompile(W_sample_up_to, (Int64, Vector{Float64}))
precompile(Wcorr_sample_up_to, (Int64, CoeffCorr{Float64}))
precompile(W_sample, (Int64, Vector{Float64}))
precompile(Wcorr_sample, (Int64, CoeffCorr{Float64}))
precompile(Wcorr_sample_mt, (Int64, Int64, CoeffCorr{Float64}))

end # module AngularCorrelation
