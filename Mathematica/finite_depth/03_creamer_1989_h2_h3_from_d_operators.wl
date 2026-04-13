(* ::Package:: *)
(*
  Re-derive Creamer et al. (1989) H2 and H3 from their D-operator notation.

  Goal
  ----
  Do not start by naming a standard Dirichlet-Neumann expansion.  Start from
  the notation used in Creamer et al. (1989):

      H = 1/2 Integral dx [ phi_s (D_z - (partial_x zeta).D_x) phi_s
                            + g zeta^2 ]                         (2.9)

  and derive the first two Hamiltonian pieces by expanding the interior
  harmonic potential around the flat surface.

  For the final Fourier check we use the 1989 deep-water symbol

      theta(k) = |k|.
*)

ClearAll["Global`*"];

Print["Creamer 1989 H2/H3 derivation from the paper's D-operator notation"];
Print[""];

Print["Paper starting point, eq. (2.9):"];
Print["  H = 1/2 Integral dx [ phi_s (D_z - (partial_x zeta).D_x) phi_s + g zeta^2 ]"];
Print[""];

Print["Let the interior harmonic potential be expanded near the flat surface:"];
Print["  g = g0 + g1 + O(zeta^2)"];
Print["  g(x,zeta(x)) = phi_s(x)"];
Print[""];

(* Boundary expansion.  The coefficient of alpha represents one power of zeta. *)
ClearAll[alpha, zeta, phiS, g0Surf, g1Surf, g0zSurf];

boundary = Normal@Series[
    g0Surf + alpha (g1Surf + zeta g0zSurf),
    {alpha, 0, 1}
  ];
boundaryOrder0 = Coefficient[boundary, alpha, 0];
boundaryOrder1 = Coefficient[boundary, alpha, 1];
g0BoundaryRule = g0Surf -> phiS;
g1BoundaryRule = Solve[boundaryOrder1 == 0, g1Surf][[1]];

Print["Boundary expansion:"];
Print["  order 0: ", boundaryOrder0 == phiS, "  ->  ", g0BoundaryRule];
Print["  order 1: ", boundaryOrder1 == 0, "  ->  ", g1BoundaryRule];
Print[""];

Print["At the flat boundary for deep water:"];
Print["  g0_z(0) = theta phi_s"];
Print["  g1(0) = -zeta theta phi_s"];
Print["  therefore g1_z(0) = theta g1(0) = -theta(zeta theta phi_s)"];
Print["  Laplace equation gives g0_zz(0) = -Delta phi_s"];
Print[""];

ClearAll[thetaPhi, thetaZetaThetaPhi, lapPhi, gradZetaDotGradPhi,
  dz0, dz1, dxCorrection, creamerOperator0, creamerOperator1,
  divZetaGradPhi];

dz0 = thetaPhi;
dz1 = -thetaZetaThetaPhi - zeta lapPhi;
dxCorrection = gradZetaDotGradPhi;

creamerOperator0 = dz0;
creamerOperator1 = FullSimplify[dz1 - dxCorrection];

divZetaGradPhi = zeta lapPhi + gradZetaDotGradPhi;

Print["Now expand Creamer's operator"];
Print["  N(zeta) phi_s = (D_z - (partial_x zeta).D_x) phi_s"];
Print["through first order:"];
Print["  N0 phi_s = ", creamerOperator0];
Print["  N1 phi_s = ", creamerOperator1];
Print["Since Div(zeta Grad phi_s) = zeta Delta phi_s + Grad[zeta].Grad[phi_s],"];
Print["  N1 phi_s = -theta(zeta theta phi_s) - Div(zeta Grad phi_s)."];
Print[""];

ClearAll[g, G0Phi, G1Phi];
hAlpha = 1/2 (
    alpha phiS (alpha G0Phi + alpha^2 G1Phi)
    + g (alpha zeta)^2
  );
h2Integrand = Coefficient[Expand[hAlpha], alpha, 2];
h3Integrand = Coefficient[Expand[hAlpha], alpha, 3];

Print["Insert N = N0 + N1 + ... into H and use eta->alpha eta, phi->alpha phi:"];
Print["  coeff[alpha^2] = ", h2Integrand];
Print["  coeff[alpha^3] = ", h3Integrand];
Print[""];

Print["Therefore:"];
Print["  H2 = 1/2 Integral dx [ g zeta^2 + phi_s theta phi_s ]"];
Print["which is Creamer et al. (1989), eq. (2.10)."];
Print[""];

Print["For H3, substitute N1 and integrate the Div term by parts:"];
Print["  H3 = 1/2 Integral dx phi_s[-theta(zeta theta phi_s) - Div(zeta Grad phi_s)]"];
Print["     = 1/2 Integral dx zeta [ |Grad phi_s|^2 - (theta phi_s)^2 ]."];
Print[""];

(* Fourier check of the H3 coefficient. *)
ClearAll[q1, q2, q3, d23, th1, th2, th3];

fromDOperators = -1/2 (d23 + th2 th3);
closureRule = d23 -> (q1 - q2 - q3)/2;
fromClosure = FullSimplify[fromDOperators /. closureRule];
deepWaterRule = {q1 -> th1^2, q2 -> th2^2, q3 -> th3^2};
fromDeepWater = FullSimplify[fromClosure /. deepWaterRule];

creamer1989H3Kernel =
  1/4 (th2^2 + th3^2 - th1^2 - 2 th2 th3);

Print["Fourier coefficient for zeta_1 phi_2 phi_3 from the derived H3:"];
Print["  before triad closure: ", fromDOperators];
Print["  after k1+k2+k3=0: ", fromClosure];
Print["  deep-water theta_i=|k_i|: ", fromDeepWater];
Print[""];
Print["Creamer et al. (1989), eq. (2.13), kernel:"];
Print["  ", creamer1989H3Kernel];
Print["Match with eq. (2.13)? ", FullSimplify[fromDeepWater == creamer1989H3Kernel]];
Print[""];

Print["Conclusion:"];
Print["  Expanding the exact operator Dz-(partial_x zeta).Dx from eq. (2.9)"];
Print["  gives H2 as eq. (2.10) and H3 as eq. (2.13)."];

