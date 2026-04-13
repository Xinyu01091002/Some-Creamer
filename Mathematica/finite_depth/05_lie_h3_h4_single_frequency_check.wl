(* ::Package:: *)
(*
  Single-frequency diagnostic for an H4/W4 correction.

  This script does not claim to have completed the full quartic normal-form
  solve.  It asks the next narrower question:

    Given the current H3-only finite-depth Creamer lambda-flow, what third-order
    cos(3 theta) surface correction would a quartic generator W4 have to supply
    in order to match the supplied finite-depth Stokes/VWA coefficient?

  The output is the target delta_C33 for the later W4 calculation.
*)

ClearAll["Global`*"];

$Assumptions = Element[sigma, Reals] && 0 < sigma < 1;

Print["H4/W4 single-frequency Stokes-limit diagnostic"];
Print[""];
Print["Convention: sigma = tanh(k h), k scaled to 1, eps = k A."];
Print[""];

creamerH3OnlyC33 =
  (3 (1 + sigma^2)*(54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 -
      4 sigma^8 - 3 sigma^10))/(64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

stokesVwaC33 =
  (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6);

residualH3Only = FullSimplify[creamerH3OnlyC33 - stokesVwaC33];
requiredW4Delta = FullSimplify[stokesVwaC33 - creamerH3OnlyC33];

ClearAll[w4DeltaC33];
solutionForCorrection = Solve[
  FullSimplify[creamerH3OnlyC33 + w4DeltaC33 - stokesVwaC33] == 0,
  w4DeltaC33
];

Print["H3-only Creamer C[3,3]:"];
Print["  ", creamerH3OnlyC33];
Print["Stokes/VWA C[3,3]:"];
Print["  ", stokesVwaC33];
Print["H3-only residual Creamer - Stokes:"];
Print["  ", residualH3Only];
Print[""];

Print["Required W4-generated surface correction delta_C[3,3]:"];
Print["  ", requiredW4Delta];
Print["Solve check:"];
Print["  ", solutionForCorrection];
Print[""];

Print["Deep-water limit sigma -> 1:"];
Print["  {Creamer_H3_only, Stokes, residual, required_delta} = ",
  FullSimplify[
    {creamerH3OnlyC33, stokesVwaC33, residualH3Only, requiredW4Delta} /. sigma -> 1
  ]
];
Print[""];

Print["Interpretation:"];
Print["  The quartic normal-form/W4 calculation must produce exactly this"];
Print["  delta_C[3,3] in the single-frequency finite-depth Stokes limit."];
Print["  This script is a target check, not the completed W4 derivation."];
