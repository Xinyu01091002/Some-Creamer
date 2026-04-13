(* ::Package:: *)
(*
  Derive the H4 source from Creamer's D-operator notation.

  This extends 03_creamer_1989_h2_h3_from_d_operators.wl by carrying the
  boundary flattening one more order:

      g = g0 + g1 + g2 + O(eta^3)
      g(x, eta(x)) = phi_s(x)
      N(eta) phi_s = (D_z - Grad[eta].D_x) phi_s
                   = N0 phi_s + N1 phi_s + N2 phi_s + ...

  The result is an operator-space H4 formula.  It is intentionally not yet a
  completed W4 normal-form calculation.
*)

ClearAll["Global`*"];

Print["Finite-depth H4 source from Creamer's D-operator notation"];
Print[""];
Print["Start from eq. (2.9) with"];
Print["  N(eta) phi_s = (D_z - Grad[eta].D_x) phi_s"];
Print["and expand the interior harmonic potential around z=0."];
Print[""];

ClearAll[eta, psi, G0, DeltaOp, gradOp, divOp];

Print["Boundary expansion g(x,eta)=psi through second order:"];
Print["  g0(0) = psi"];
Print["  g1(0) = -eta G0 psi"];
Print["  g2(0) = -eta g1_z(0) - 1/2 eta^2 g0_zz(0)"];
Print["        =  eta G0(eta G0 psi) + 1/2 eta^2 Delta psi"];
Print[""];

Print["Flat-boundary derivative rules:"];
Print["  g0_z   = G0 psi"];
Print["  g0_zz  = -Delta psi"];
Print["  g0_zzz = -Delta G0 psi"];
Print["  g1_z   = -G0(eta G0 psi)"];
Print["  g1_zz  =  Delta(eta G0 psi)"];
Print["  g2_z   =  G0[eta G0(eta G0 psi) + 1/2 eta^2 Delta psi]"];
Print[""];

Print["Expand"];
Print["  D_z psi = g_z(x,eta)"];
Print["          = g0_z + g1_z + g2_z + eta(g0_zz+g1_zz) + 1/2 eta^2 g0_zzz + ..."];
Print["  D_x psi = g_x(x,eta)"];
Print["          = g0_x + g1_x + eta g0_xz + ..."];
Print["so N = D_z - Grad[eta].D_x."];
Print[""];

Print["The resulting operator pieces are:"];
Print["  N0 psi = G0 psi"];
Print["  N1 psi = -G0(eta G0 psi) - Div(eta Grad psi)"];
Print[""];

Print["  N2 psi ="];
Print["    G0[eta G0(eta G0 psi)]"];
Print["    + 1/2 G0[eta^2 Delta psi]"];
Print["    + eta Delta(eta G0 psi)"];
Print["    - 1/2 eta^2 Delta(G0 psi)"];
Print["    + |Grad eta|^2 G0 psi"];
Print[""];
Print["Equivalently, using"];
Print["  1/2 Delta(eta^2 G0 psi)"];
Print["    = eta Delta(eta G0 psi) - 1/2 eta^2 Delta(G0 psi) + |Grad eta|^2 G0 psi,"];
Print["we can write the compact form"];
Print["  N2 psi = G0 eta G0 eta G0 psi"];
Print["           + 1/2 G0(eta^2 Delta psi)"];
Print["           + 1/2 Delta(eta^2 G0 psi)."];
Print[""];

ClearAll[alpha, phiS, zeta, g, G0Phi, G1Phi, G2Phi];
hAlpha = 1/2 (
    alpha phiS (alpha G0Phi + alpha^2 G1Phi + alpha^3 G2Phi)
    + g (alpha zeta)^2
  );

h2Integrand = Coefficient[Expand[hAlpha], alpha, 2];
h3Integrand = Coefficient[Expand[hAlpha], alpha, 3];
h4Integrand = Coefficient[Expand[hAlpha], alpha, 4];

Print["Alpha bookkeeping in H(alpha eta, alpha phi_s):"];
Print["  coeff[alpha^2] = ", h2Integrand];
Print["  coeff[alpha^3] = ", h3Integrand];
Print["  coeff[alpha^4] = ", h4Integrand];
Print[""];
Print["Therefore:"];
Print["  H4 = 1/2 Integral dx phi_s N2(eta,eta) phi_s"];
Print["with N2 as written above."];
Print[""];

Print["Finite-depth Fourier rule for later use:"];
Print["  G0(k) = theta(k) = |k| tanh(|k| h)"];
Print["  Delta acting on mode k gives -|k|^2"];
Print[""];

ClearAll[th3, th4, th23, q3, q4];
unsymH4Kernel =
  1/2 th4 th23 th3 - 1/4 (th4 q3 + q4 th3);

Print["For one unsymmetrized monomial eta_1 eta_2 phi_3 phi_4,"];
Print["with k1+k2+k3+k4=0 and k23=k2+k3, the compact N2 form gives"];
Print["  H4 unsym kernel = ", unsymH4Kernel];
Print["where th_i=theta(k_i), th23=theta(k2+k3), and q_i=|k_i|^2."];
Print["The symmetrized quartic kernel should average over the eta and phi slots"];
Print["before it is used in a W4 normal-form solve."];
