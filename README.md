# AngularCorrelation.jl

A package to calculate angular distributions and correlations for γ cascades.
The direction and polarization of the zeroth γ-ray are known.
It is aligned with the z axis, the xz plane is the plane of (linear) polarization.

The angular correlation of subsequently emitted photons can be obtained
- … for any single emitted photon after (optionally) an arbitrary long cascade of unobserved intermediate photons.
  This is also referred to as *angular distribution*.
- … between the first emitted photon and any (single) subsequent photon after (optionally) an arbitrary long cascade of unobserved intermediate photons.

## Installation

In Julia, run

```julia
import Pkg; Pkg.add(url="https://github.com/op3/AngularCorrelation.jl.git")
```

## Usage

The syntax is straightforward,
structs are used to provide an abstract representation of the states and transitions participating in the cascade.

The following code can be used to calculate a 0⁺ → 1⁺ → 2 → 0 cascade:

```julia
using AngularCorrelation

S0 = State(0, Parity.positive)
g0 = M1()
S1 = State(1, Parity.positive)
g1 = Dipole()
S2 = State(2, Parity.unknown)
g2 = Quadrupole()
S3 = State(0, Parity.unknown)

theta1 = π/2
phi1 = π/2
theta2 = π/2
phi2 = 3π/2

Wcorr(theta1, phi1, theta2, phi2,
      S0, g0, S1, g1, S2, g2, S3)
```

Supplying parity and electromagnetic character information for anything
but the first transition and first two states is useless.
The parity of the states is optional (but a useful cross-check to make sure the transition makes sense)

Angular distributions are calculated in the following way:

```julia
theta = π/2
phi = π/2
W(theta, phi, S0, g0, S1, g1, S2, g2, S3)
```

The cascade can be arbitrarily long.
The function `W` returns the angular distribution of the last emitted photon in the given cascade.

For states with half-integer angular momentum (i.e., odd-mass nuclei),
the angular momentum has to be given as a Rational number:

```julia
s32 = State(3//2)
```

Multipole mixing ratios δ are optional arguments of the transitions:

```julia
mixed = M1(-0.34)
```

## Conventions

The Integral over the complete probability distribution
( ∫∫∫∫ Wcorr(θ₁, ϕ₁,θ₂, ϕ₂) sin(θ₁) sin(θ₂) dθ₁ dφ₁ dθ₂ dφ₂ and
∫∫ W(θ, ϕ) sin(θ) dθ dφ )
is equal to 4π.
Thus, if the direction of the first or second γ is fixed,
one obtains an angular distribution for the other γ.

The KSW convention by Krane, Steffen, Wheeler [\[2\]](#ref-2) for the multipole mixing ratio δ is used.
It is assumed that only the first and second-order multipolarities contribute.

## License<a name="license"></a>

© O. Papst [`<opapst@ikp.tu-darmstadt.de>`](mailto:opapst@ikp.tu-darmstadt.de)

AngularCorrelation.jl is distributed under the terms of the GNU General Public License, version 3 or later.
See the [`LICENSE`](LICENSE) file.

## References

<a name="ref-1">[1]</a> R. M. Steffen and K. Alder, “Angular distribution and correlation of gamma rays”, in *The electromagnetic interaction in nuclear spectroscopy*, edited by W. D. Hamilton (North-Holland, Amsterdam, 1975) Chap. 12, pp. 505–582, ISBN: 978-0-4441-0519-6.

<a name="ref-2">[2]</a> K. S. Krane, R. M. Steffen and R. M. Wheeler, “Directional correlations of gamma radiations emitted from nuclear states oriented by nuclear reactions or cryogenic methods”, At. Data Nucl. Data Tables **11**, 351 (1973). [`doi:10.1016/S0092-640X(73)80016-6`](https://doi.org/10.1016/S0092-640X(73)80016-6).  

## Acknowledgments

The author would like to thank

- J. Isaak for providing a reference implementation to test my results with
- U. Friman-Gayer for valuable discussion. In particular, the author would like to advertise his [alpaca](https://github.com/u-eff-gee/alpaca) code
  that can be used to calculate angular correlations between two γ-rays
  and also sample cascades of arbitrary length.
  His code uses the Biedenharn convention, though.

This work has been funded by the State of Hesse under the grant “Nuclear Photonics” within the LOEWE program.
