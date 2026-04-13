(* ::Package:: *)
(*
  Check the H2/H3 Hamiltonian expansion used before attempting H4.

  This script starts from the Hamiltonian in Creamer et al. (1989), eq. (2.9),
  written with the Dirichlet-Neumann operator

      H = 1/2 Integral [ g eta^2 + phi G(eta) phi ] dx
      G(eta) = G0 + G1(eta) + ...

  with finite-depth flat-bottom symbol

      G0(k) = theta(k) = |k| Tanh[|k| h].

  The first-order DNO correction is

      G1(eta) phi = -G0 eta G0 phi - Div[eta Grad[phi]].

  Therefore

      H2 = 1/2 Integral [ g eta^2 + phi G0 phi ] dx
      H3 = 1/2 Integral eta [ |Grad phi|^2 - (G0 phi)^2 ] dx.

  In Fourier variables and under triad closure k1+k2+k3=0, this should
  reproduce Wright & Creamer (1994), eq. (10):

      H3 kernel = 1/4 ( |k2|^2 + |k3|^2 - |k1|^2 - 2 theta2 theta3 ).

  The deep-water limit theta_i -> |k_i| then reproduces Creamer et al. (1989),
  eq. (2.13).
*)

ClearAll["Global`*"];

Print["Hamiltonian H2/H3 check from the finite-depth DNO expansion"];
Print[""];

Print["Start from Creamer et al. (1989), eq. (2.9):"];
Print["  H(eta,phi) = 1/2 Integral dx [ phi (Dz - Grad[eta].Dx) phi + g eta^2 ]"];
Print[""];
Print["Define the Dirichlet-Neumann operator:"];
Print["  G(eta) = Dz - Grad[eta].Dx"];
Print["This is the same non-unit-normal DNO because Creamer's D_j are the"];
Print["interior derivatives evaluated at z=eta, and eq. (2.7) gives"];
Print["  Dx + Grad[eta] Dz = Grad"];
Print["so Dz - Grad[eta].Dx is the usual phi_z - Grad[eta].Grad_x phi"];
Print["written in surface variables."];
Print["so the same Hamiltonian is:"];
Print["  H(eta,phi) = 1/2 Integral dx [ phi G(eta) phi + g eta^2 ]"];
Print[""];

Print["Finite-depth DNO expansion used:"];
Print["  G0(k) = theta(k) = |k| Tanh[|k| h]"];
Print["  G1(eta) phi = -G0 eta G0 phi - Div[eta Grad[phi]]"];
Print["One way to see G1 is to flatten the boundary at z=0:"];
Print["  phi = phi0 + phi1 + ..."];
Print["  phi0(0)=surface phi,     phi0_z(0)=G0 phi"];
Print["  phi1(0)=-eta G0 phi,     phi1_z(0)=-G0(eta G0 phi)"];
Print["  evaluating Dz at z=eta adds eta phi0_zz(0)=-eta Laplacian phi"];
Print["  subtracting Grad[eta].Dx adds -Grad[eta].Grad[phi]"];
Print["  so -eta Laplacian phi - Grad[eta].Grad[phi] = -Div[eta Grad[phi]]."];
Print[""];

Print["Introduce a bookkeeping parameter alpha:"];
Print["  eta -> alpha eta"];
Print["  phi -> alpha phi"];
Print["  G(alpha eta) = G0 + alpha G1(eta) + O(alpha^2)"];
Print[""];

(* Symbolic noncommutative placeholders for the operator-level truncation.
   This is not yet a Fourier calculation; it only tracks variable degree. *)
ClearAll[alpha, eta, phi, g, G0Phi, G1Phi, gradPhi2, g0Phi2];

hAlpha = 1/2 (alpha phi (G0Phi alpha + alpha^2 G1Phi) + g (alpha eta)^2);
h2Alpha = Coefficient[Expand[hAlpha], alpha, 2];
h3Alpha = Coefficient[Expand[hAlpha], alpha, 3];

Print["Alpha-truncated integrand before using the explicit form of G1:"];
Print["  coeff[alpha^2] = ", h2Alpha];
Print["  coeff[alpha^3] = ", h3Alpha];
Print[""];

Print["Interpreting those coefficients gives:"];
Print["  H2 = 1/2 Integral dx [ g eta^2 + phi G0 phi ]"];
Print["  H3 = 1/2 Integral dx phi G1(eta) phi"];
Print[""];

Print["Using G1 and integrating the Div term by parts:"];
Print["  Integral phi[-Div(eta Grad phi)] dx = Integral eta |Grad phi|^2 dx"];
Print["  Integral phi[-G0 eta G0 phi] dx = -Integral eta (G0 phi)^2 dx"];
Print["so:"];
Print["  H3 = 1/2 Integral dx eta [ |Grad phi|^2 - (G0 phi)^2 ]"];
Print[""];

(* Let q_i = |k_i|^2, theta_i = theta(k_i), and d23 = k2 dot k3. *)
ClearAll[q1, q2, q3, d23, th1, th2, th3];

fromDNO = -1/2 (d23 + th2 th3);

closureRule = d23 -> (q1 - q2 - q3)/2;
fromDNOClosure = FullSimplify[fromDNO /. closureRule];

wrightCreamer1994 = 1/4 (q2 + q3 - q1 - 2 th2 th3);

Print["Fourier H3 coefficient from H3 = 1/2 eta[(grad phi)^2 - (G0 phi)^2]:"];
Print["  before triad closure: ", fromDNO];
Print["  using k1+k2+k3=0, so k2.k3=(|k1|^2-|k2|^2-|k3|^2)/2:"];
Print["  from DNO:              ", fromDNOClosure];
Print["  Wright-Creamer (1994): ", wrightCreamer1994];
Print["  Match? ", FullSimplify[fromDNOClosure == wrightCreamer1994]];
Print[""];

deepWaterRules = {q1 -> th1^2, q2 -> th2^2, q3 -> th3^2};
deepWaterFromDNO = FullSimplify[fromDNOClosure /. deepWaterRules];
creamer1989 = 1/4 (th2^2 + th3^2 - th1^2 - 2 th2 th3);

Print["Deep-water specialization |k_i| = theta_i:"];
Print["  from finite-depth formula: ", deepWaterFromDNO];
Print["  Creamer et al. (1989):     ", creamer1989];
Print["  Match? ", FullSimplify[deepWaterFromDNO == creamer1989]];
Print[""];

Print["Takeaway:"];
Print["  The finite-depth H2/H3 pieces are reproduced by the standard DNO"];
Print["  expansion through G1. This is the low-order check to pass before"];
Print["  deriving or testing any H4/W4 correction."];
