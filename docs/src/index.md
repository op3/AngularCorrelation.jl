# AngularCorrelation.jl

A package to calculate angular distributions and correlations for γ cascades.
The direction and polarization of the zeroth γ-ray are known.
It is aligned with the z axis, the xz plane is the plane of (linear) polarization.

The angular correlation of subsequently emitted photons can be obtained
- … for any single emitted photon after (optionally) an arbitrary long cascade of unobserved intermediate photons.
  This is also referred to as *angular distribution*.
- … between the first emitted photon and any (single) subsequent photon after (optionally) an arbitrary long cascade of unobserved intermediate photons.

## Definition of cascades

```@docs
Transition
State
Parity
EMCharacter
AngularCorrelation.alt_char
E1
M1
E2
M3
Dipole
Quadrupole
```

## Angular distributions

```@docs
W_coeff
W
W_sample
W_sample_up_to
```

## Angular correlations


```@docs
Wcorr_coeff
Wcorr
Wcorr_sample
Wcorr_sample_mt
Wcorr_sample_up_to
```

# Coefficients

```@docs
AngularCorrelation.F
AngularCorrelation.A
AngularCorrelation.B
AngularCorrelation.B_lpol
AngularCorrelation.U
AngularCorrelation.kappa
```
