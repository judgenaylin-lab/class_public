# TQS scalar second variation on a moving homogeneous mesh background

This note supersedes the pinned-background scalar closure when dot(phi_bar) != 0. It is derived within the current production truncation: the standard CLASS LambdaCDM expansion is retained as an external background while the homogeneous TQS mesh is evolved self-consistently on it. The full TQS Friedmann normalization is a separate extension.

## Definitions

Use grid cosmic time t with v = dot(phi_bar), F = A Xi, and D = k^2/a^2.
Let u = delta phi_U be the unitary-clock mesh perturbation and w1 = delta(D_T phi) = dot(u) - v n_hat.

The homogeneous force is

F_TQS = 3 beta_phi Gamma_b H^2 + 3 beta_phi sqrt(A) (rho+3p) + m_L^2 s_C phi_bar.

The homogeneous mesh equation is

dot(v) + (3H + dot(Xi)/Xi) v - 2 beta_phi v^2 + F_TQS/(A Xi) = 0.

## Quadratic kinetic contribution

For S_kin = (M_*^2/2) integral N sqrt(gamma) F(phi) (D_T phi)^2, with F_phi/F = -4 beta_phi, the moving-background quadratic kinetic sector produces the compact perturbations

delta rho_mesh/M_*^2 = F v w1 - 2 beta_phi F v^2 u + m_L^2 s_C phi_bar u,

and a mesh momentum contribution F v u to the scalar 0i constraint.

## Corrected unitary constraints

2 Gamma q_hat + B s_hat + 2 Gamma_phi H u + F v u = -J_hat,

4 D zeta_hat + 4 Gamma H (3 q_hat+s_hat) + 6 Gamma_phi H^2 u + 2 alpha_K D n_hat - 2 delta rho_mesh/M_*^2 = 2 R_hat.

These reduce to the pinned equations when v -> 0.

## Physical Newtonian dictionary

CLASS psi = tilde Phi and CLASS phi = tilde Psi.

The physical Newtonian mesh perturbation is delta phi_N = u + v pi.

n_hat = tilde Phi - beta_phi (u+v pi) - dot(pi),
zeta_hat = -tilde Psi + beta_phi (u+v pi) - H pi,
s_hat = D pi/A.

With q_hat = -X,

X = J_hat/(2 Gamma) + B D pi/(2 Gamma A) + (Gamma_phi/Gamma) H u + (F v/(2 Gamma)) u.

The physical metric evolution is

dot(tilde Psi) = -H tilde Phi - dot(H) pi + X + beta_phi [ w1 + v tilde Phi + (H-beta_phi v)u + (dot(v)+H v-beta_phi v^2)pi ].

The Hamiltonian constraint determines n_hat algebraically and dot(pi) = tilde Phi - beta_phi(u+v pi) - n_hat.

The slip remains GR-form because A t_hat = D pi even when A=A(t), while the homogeneous mesh has no linear anisotropic stress.

## First-order moving-background mesh system

Define M_move^2 = partial_phi(F_TQS) + 4 beta_phi F_TQS, with

partial_phi(F_TQS) = m_L^2 s_C + 3 beta_phi^2 Gamma_bb H^2 - 6 beta_phi^2 sqrt(A) (rho+3p).

The kinematic equation is dot(u) = w1 + v n_hat.

The velocity perturbation obeys

dot(w1) = -(3H + dot(Xi)/Xi - 4 beta_phi v) w1 - v deltaK + (dot(v)+dot(Xi)v/Xi) n_hat - (mu D + M_move^2)u/F - [2 beta_phi Gamma_b H deltaK + 3 beta_phi sqrt(A)(delta rho+3 delta p)]/F,

with deltaK = -3X + D pi/A.

Using the homogeneous equation, dot(v)+dot(Xi)v/Xi = -3Hv + 2 beta_phi v^2 - F_TQS/F, so no numerical second derivative of the background field is needed.

## CLASS state convention

The existing scalar slots are retained:
- tqs_varphi = u = delta phi_U
- tqs_varphi_prime = a w1, not u' on a moving background.

This avoids n_hat' and pi'' while reducing continuously to the pinned implementation when phi_bar' -> 0.

## Scope

This is the scalar second variation of the currently implemented matched-clock + quadratic mesh-kinetic + regular-pinning truncation on a moving homogeneous mesh background, with s_C kept as the homogeneous saturated FLRW response. Perturbations of s_C and the full TQS Friedmann equation are not included in this stage.