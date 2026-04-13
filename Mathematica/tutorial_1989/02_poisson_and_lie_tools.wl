(* Creamer 1989 tutorial, section 2 -> 3 bridge.
   Goal:
     Explain how H2 generates the linear equations and why the Poisson bracket
     matters before the Lie transform is introduced.
*)

ClearAll["Global`*"];

$Assumptions = {
  g > 0,
  k > 0,
  Element[{g, k}, Reals]
};

thetaDeep[k_] := Abs[k];
omegaDeep[k_] := Sqrt[g Abs[k]];

Print["Section 2 -> 3 bridge: canonical equations and Poisson bracket"];
Print[""];
Print["Paper-level Poisson bracket, equation (3.6):"];
Print["  {A,B} = Integral d^2x"];
Print["           [ dA/dzeta(x) dB/dphi_s(x) - dA/dphi_s(x) dB/dzeta(x) ]"];
Print[""];
Print["Why is this useful?"];
Print["  Because zeta and phi_s are canonical variables, so once H is known, the"];
Print["  equations of motion follow immediately from the Hamiltonian structure."];
Print[""];
Print["Canonical equations:"];
Print["  zeta_t   =   dH/dphi_s"];
Print["  phi_s,t  = - dH/dzeta"];
Print[""];
Print["If only H2 is kept, then in Fourier space one gets"];
Print["  zetaDot_k = |k| phi_k"];
Print["  phiDot_k  = - g zeta_k"];
Print[""];
Print["Eliminating phi_k gives"];
Print["  zetaDDot_k + g |k| zeta_k = 0"];
Print["which is the deep-water linear dispersion relation"];
Print["  omega_k^2 = g |k|"];
Print[""];
Print["So H2 is called the quadratic or linear Hamiltonian because by itself it"];
Print["generates the standard linear gravity-wave dynamics."];
Print[""];
Print["Why does section 3 introduce a Lie transform?"];
Print["  Because the paper wants a canonical change of variables that removes H3"];
Print["  while preserving the Hamiltonian structure."];
Print[""];
Print["Equation (3.7) then introduces an auxiliary parameter lambda through"];
Print["  partial_lambda A = {A, W}"];
Print["where A is either Z or Phi."];
Print["Lambda is not physical time. It parameterizes the canonical deformation from"];
Print["the original variables to the transformed variables."];
