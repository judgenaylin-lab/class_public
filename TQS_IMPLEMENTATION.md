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
