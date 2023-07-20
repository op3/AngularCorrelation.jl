# AngularCorrelation.jl

A package to calculate angular correlations for triple γ cascades.
The direction and polarization of the zeroth γ are known.
It is aligned with the z axis, the xz plane is the plane of (linear) polarization.

This is code can be used to calculate the angular correlation of
the first two γ-rays emitted from a nucleus that was excited
by a linearly polarized γ-ray beam (the first γ).

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

Supplying polarity and electromagnetic character information for anything
but the first transition and first two states is useless.

## Conventions

The Integral over the complete probability distribution ( ∫∫∫∫ Wcorr(θ₁, ϕ₁,θ₂, ϕ₂) sin(θ₁) sin(θ₂) dθ₁ dφ₁ dθ₂ dφ₂ )
is equal to 4π.
Thus, if the direction of the first or second γ is fixed,
one obtains an angular distribution for the other γ.

The KSW convention by Krane, Steffen, Wheeler [\[2\]](#ref-2) for the multipole mixing ratio δ is used.

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
