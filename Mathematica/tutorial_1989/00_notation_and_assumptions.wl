(* Creamer 1989 tutorial, section 0.
   Purpose:
     Fix the paper's notation before any derivation begins.

   This file is a dictionary, not yet a derivation.
*)

ClearAll["Global`*"];

$Assumptions = {
  g > 0,
  k > 0,
  q > 0,
  Element[{g, h, k, q, x, y, t}, Reals]
};

thetaDeep[k_] := Abs[k];
omegaDeep[k_] := Sqrt[g Abs[k]];
hilbertSymbol1D[k_] := -I Sign[k];
thetaFinite[k_, h_] := k Tanh[k h];

Print["Section 0: notation"];
Print[""];
Print["Paper-level canonical variables:"];
Print["  zeta(x,t)   = free-surface elevation"];
Print["  phi_s(x,t)  = velocity potential evaluated at the free surface"];
Print[""];
Print["Deep-water operator from equation (2.11):"];
Print["  theta(k) = ", thetaDeep[k]];
Print[""];
Print["Deep-water linear dispersion relation:"];
Print["  omega_k = ", omegaDeep[k], " and omega_k^2 = g |k|"];
Print["  Later this will come directly from the quadratic Hamiltonian H2."];
Print[""];
Print["1D Hilbert-transform convention used in section 4:"];
Print["  H[f]_k = ", hilbertSymbol1D[k], " f_k"];
Print["  For q > 0, H[Cos(q x)] = ", Sin[q x]];
Print[""];
Print["Finite-depth comparison symbol, only for later reference:"];
Print["  theta_h(k) = ", thetaFinite[k, h]];
