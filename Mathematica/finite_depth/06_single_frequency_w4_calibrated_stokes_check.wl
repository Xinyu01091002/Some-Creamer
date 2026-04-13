(* ::Package:: *)
(*
  Calibrated single-frequency W4 correction from the finite-depth Stokes limit.

  This script implements the narrow check requested after the H4 source was
  identified:

    1. Use the current H3-only finite-depth Creamer single-frequency C[3,3].
    2. Use the supplied finite-depth Stokes/VWA third-harmonic coefficient.
    3. Define the one-mode W4 surface correction required to make them match.
    4. Verify the corrected coefficient equals Stokes and vanishes back to the
       deep-water correction in the sigma -> 1 limit.

  Important: this is a calibrated one-mode target for W4.  It does not yet
  derive the complete quartic normal-form generator W4(k1,k2,k3,k4) from
  K4 = H4 + 1/2 {H3,W3}.  It is the acceptance test that such a W4 derivation
  must satisfy in the finite-depth regular-wave limit.
*)

ClearAll["Global`*"];

$Assumptions = Element[sigma, Reals] && 0 < sigma < 1;

Print["Calibrated single-frequency W4 correction from Stokes/VWA C[3,3]"];
Print[""];
Print["Convention:"];
Print["  sigma = tanh(k h), k scaled to 1, eps = k A"];
Print["  eta(theta) contains C[3,3] eps^3 Cos[3 theta] at third order."];
Print[""];

(* H3-only result from 01_single_frequency_limit.wl after T[n] -> n tanh(n atanh(sigma)). *)
creamerH3OnlyC33 =
  (3 (1 + sigma^2) (54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 -
      4 sigma^8 - 3 sigma^10))/
    (64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

stokesVwaC33 =
  (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6);

h3OnlyResidual = FullSimplify[creamerH3OnlyC33 - stokesVwaC33];

(* The calibrated W4 contribution to the coordinate transform must be the
   negative of the H3-only residual in this single-frequency observable. *)
w4CalibratedDeltaC33 = FullSimplify[-h3OnlyResidual];
correctedC33 = FullSimplify[creamerH3OnlyC33 + w4CalibratedDeltaC33];
correctedResidual = FullSimplify[correctedC33 - stokesVwaC33];

(* Minimal one-mode generator ansatz.
   Let z_1 = eps/2.  A generator term W4_min = w z_1^3 phi_-3 + c.c.
   gives a complex z_3 correction proportional to w eps^3/8.  Since the
   real cos(3 theta) coefficient is twice the positive complex coefficient,
   delta_C33 = w/4. *)
w4OneModeGeneratorCoeff = FullSimplify[4 w4CalibratedDeltaC33];

Print["H3-only Creamer C[3,3]:"];
Print["  ", creamerH3OnlyC33];
Print["Stokes/VWA C[3,3]:"];
Print["  ", stokesVwaC33];
Print["H3-only residual:"];
Print["  ", h3OnlyResidual];
Print[""];

Print["Calibrated one-mode W4 delta_C[3,3]:"];
Print["  ", w4CalibratedDeltaC33];
Print["Minimal one-mode generator coefficient w in W4_min = w z_1^3 phi_-3 + c.c.:"];
Print["  ", w4OneModeGeneratorCoeff];
Print["Corrected C[3,3] = H3-only + W4 delta:"];
Print["  ", correctedC33];
Print["Corrected residual vs Stokes:"];
Print["  ", correctedResidual];
Print[""];

Print["Deep-water limit sigma -> 1:"];
Print["  {H3-only, W4 delta, corrected, Stokes, corrected residual} = ",
  FullSimplify[
    {creamerH3OnlyC33, w4CalibratedDeltaC33, correctedC33,
      stokesVwaC33, correctedResidual} /. sigma -> 1
  ]
];
Print[""];

Print["Finite-depth example sigma=tanh(2):"];
Print["  {H3-only, W4 delta, corrected, Stokes} = ",
  N[
    {creamerH3OnlyC33, w4CalibratedDeltaC33, correctedC33,
      stokesVwaC33} /. sigma -> Tanh[2],
    14
  ]
];
Print[""];

Print["Acceptance target for a future full W4 derivation:"];
Print["  The W4 transform contribution to eta_3 in this single-frequency limit"];
Print["  must reduce to the calibrated delta_C[3,3] above."];
