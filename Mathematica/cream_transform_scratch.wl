(* Minimal Wolfram Language scratchpad for Creamer-transform derivations. *)

ClearAll["Global`*"];

$Assumptions = {
  A > 0,
  k > 0,
  eps > 0,
  Element[{x, t, h}, Reals]
};

theta[x_, t_] := k x - w t;

(* Toy model: Eulerian reconstruction from horizontal remapping. *)
xi[x_, t_] := eps Sin[theta[x, t]];
etaEuler[x_, t_] := A Cos[theta[x - xi[x, t], t]];

Print["Series for remapped profile through O(eps^2):"];
Print[Normal @ Series[etaEuler[x, t], {eps, 0, 2}] // TrigReduce // Simplify];

(* Finite-depth linear operator symbol often appearing in the literature. *)
thetaFD[k_, h_] := k Tanh[k h];

Print["Finite-depth operator symbol thetaFD(k,h):"];
Print[thetaFD[k, h]];

Print["Small-kh expansion:"];
Print[Normal @ Series[thetaFD[k, h], {h, 0, 5}] // Simplify];
