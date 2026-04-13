(* Scratchpad for temporary experiments.

   Rule for this file:
     Even though this is a scratch file, do not jump directly to a simplified formula
     without labeling what was simplified.

   Minimal pattern to follow:
     1. Write the original paper-level object or quote the target equation.
     2. State the simplification being made:
        - one mode
        - periodic interval
        - small parameter expansion
        - toy canonical pair
        - finite-harmonic truncation
     3. State why the simplification is useful.

   This file is intentionally lightweight, but it should still stay honest about the
   bridge from the paper to the local experiment.
*)

ClearAll["Global`*"];

$Assumptions = {
  A > 0,
  k > 0,
  eps > 0,
  Element[{x, t, h}, Reals]
};

(* Example scratch experiment.
   Original geometric idea from the 1D discussion:
     nonlinear geometry can be viewed as a horizontal remapping.
   Simplification used here:
     choose an explicit small horizontal shift xi = eps Sin(kx-wt).
   Reason:
     quickly inspect the harmonic content created by the remapping. *)

theta[x_, t_] := k x - w t;
xi[x_, t_] := eps Sin[theta[x, t]];
etaEuler[x_, t_] := A Cos[theta[x - xi[x, t], t]];

Print["Series for remapped profile through O(eps^2):"];
Print[Normal @ Series[etaEuler[x, t], {eps, 0, 2}] // TrigReduce // Simplify];

(* Finite-depth symbol kept here for quick comparisons. *)
thetaFD[k_, h_] := k Tanh[k h];

Print["Finite-depth operator symbol thetaFD(k,h):"];
Print[thetaFD[k, h]];

Print["Small-kh expansion:"];
Print[Normal @ Series[thetaFD[k, h], {h, 0, 5}] // Simplify];
