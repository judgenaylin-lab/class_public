# TQS production implementation notes

This branch integrates the current regular-pinned, matched-clock TQS scalar
sector directly into CLASS's Newtonian-gauge perturbation module.

## Variable map

CLASS uses

```
ds^2 = a^2[-(1+2 psi)d tau^2 + (1-2 phi) dx^2]
```

and the TQS physical Newtonian variables are mapped as

```
CLASS psi = tilde Phi
CLASS phi = tilde Psi
```

The TQS scalar state adds

```
pi             clock Stückelberg displacement
varphi         mesh perturbation
varphi_prime   d varphi / d tau
```

The physical slip relation remains

```
k^2 (tilde Phi - tilde Psi)
 = 12 pi G_* a_phys^2 sum_s (rho_s+p_s) sigma_s
```

which has the same CLASS form as GR on the matched pinned branch.

## Current production assumptions

- LambdaCDM background expansion is retained for the first CLASS likelihood pass.
- A is homogeneous and constant.
- s_C follows the saturated FLRW equilibrium response
  H^2/(H^2+H_C^2).
- regular pinning uses m_eff^2=m_L^2 s_C.
- photons, baryons, CDM and neutrino hierarchies are native CLASS.
- TQS modifies only the scalar metric/clock/mesh closure.

## Validation order

1. compile CLASS;
2. reproduce baseline LambdaCDM with tqs_enable=no;
3. run TQS without NaN/ODE failure;
4. verify beta_phi -> 0 and A -> 1 limits;
5. compare TT/TE/EE/phiphi;
6. increase precision and l_max;
7. only then attach observational likelihoods.

Do not interpret the smoke-test spectra as a fit before these checks pass.


## Dynamic-background production result (2026-10-01)

The fixed-background implementation was found to be inconsistent for nonzero
beta_phi.  The homogeneous mesh is now evolved as

```
phi_TQS'' + (2 Hconf + Xi'/Xi) phi_TQS'
           - 2 beta_phi (phi_TQS')^2
           + a^2 F_TQS/(A Xi) = 0
```

with

```
F_TQS =
  3 beta_phi Gamma_b H^2
  + 3 beta_phi sqrt(A) (rho + 3 p)
  + m_L^2 s_C phi_TQS
```

and

```
A(phi_TQS) = A_ref exp(-4 beta_phi phi_TQS).
```

This first stage intentionally keeps the standard CLASS LambdaCDM Friedmann
background so that homogeneous mesh consistency can be tested separately from
the not-yet-final TQS Friedmann equation.

For A_ref=0.98, beta_phi=0.111272 and m_L/H0=1e5:

- mu_C <= 1800 leaves the healthy branch A<1.
- mu_C = 2000 survives but reaches A_max ~= 0.999257.
- mu_C = 2200 gives A_max ~= 0.997612.
- mu_C = 2500 gives A_rec ~= 0.995361 and A_max ~= 0.995614.
- mu_C = 3000 gives A_rec ~= 0.992902 and A_max ~= 0.993136.

The current provisional production smoke point is mu_C=2500.

At that point the native CLASS CMB integration completes.  The first
background-fed perturbation comparison gives percent-level spectral shifts
(roughly +4% TT/EE and +8% lensing over much of the acoustic range).

IMPORTANT: these are diagnostic spectra, not final TQS predictions.  Once
phi_TQS(t) is non-constant, the exact second variation contains additional
terms proportional to the homogeneous phi_TQS' and phi_TQS''.  Those terms
have not yet been derived/implemented.  A precision likelihood must wait for
that moving-background second variation.
