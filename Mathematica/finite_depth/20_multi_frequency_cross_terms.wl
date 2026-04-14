(* ::Package:: *)
(*  Multi-frequency cross-interaction check: H3-only vs H3+H4.

    Strategy
    --------
    PART A  uses modes = {-3,...,-1, 1,...,3} to compute cross terms whose
            output mode lies within the standard truncation, e.g.
              (-1, 2, 2) -> 3   eps1 eps2^2 at output mode 3
              (-1, 1, 3) -> 3   eps1^2 eps3 (trichromatic, needs mode 3 active)
              (-2, 2, 3) -> 3   eps2^2 eps3 (trichromatic)

    PART B  extends modes to {-6,...,-1, 1,...,6} so that output modes 4, 5, 6
            become reachable.  W3, H4, W4Total are rebuilt on this larger set.
            Cross terms covered:
              (1,  1, 2) -> 4   eps1^2 eps2
              (1,  2, 2) -> 5   eps1 eps2^2
              (2,  2, 2) -> 6   eps2^3  (single-freq check for mode 2)
              (1,  1, 3) -> 5   eps1^2 eps3
              (1,  2, 3) -> 6   eps1 eps2 eps3
              (2,  2, 3) -> 7   (not in new range, left as note)

    Coefficient convention: 2*C so that a real wave of amplitude A gives
    eta_nl[m] ~ C A^3 cos(m x), matching the Stokes C33 normalisation.

    Variable names: camelCase only; no underscores (Mathematica reserved).
*)

ClearAll["Global`*"];

$Assumptions = (And @@ Table[T[n] > 0, {n, 1, 12}]) &&
  Element[sigma, Reals] && 0 < sigma < 1;

(* ------------------------------------------------------------------ *)
(*   Kernel helpers (independent of mode set)                          *)
(* ------------------------------------------------------------------ *)

denBase[t1_, t2_, t3_] :=
  t1 (t2 + t3 - t1) + t2 (t1 + t3 - t2) + t3 (t1 + t2 - t3);

dFinite[t1_, t2_, t3_, q1_, q2_, q3_] := Module[{den, num},
  den = denBase[t1, t2, t3];
  If[t1 === 0 || t2 === 0 || t3 === 0 || TrueQ[den == 0], 0,
    num = q1 (t1^2 - (t2-t3)^2) + q2 (t2^2-(t1-t3)^2) + q3 (t3^2-(t1-t2)^2)
        - 2 t1 t2 t3 (t1 + t2 + t3);
    num/(12 den)]];

bFinite[t1_, t2_, t3_, q1_, q2_, q3_] := Module[{den, num},
  den = denBase[t1, t2, t3];
  If[TrueQ[den == 0], 0,
    num = t3 (t3 (t1+t2) - t1^2 - t2^2)
        + t3 (q1+q2-2 q3) + (t1-t2)(q1-q2);
    num/(2 den)]];

thetaOf[n_Integer] /; n != 0 := T[Abs[n]];
thetaOf[0] := 0;
omegaOf[n_Integer] /; n != 0 := Sqrt[T[Abs[n]]];
ksqOf[n_Integer] := n^2;

(* cached kernels: use global cache, keyed by mode triple *)
DkerVal[a_, b_, c_] := DkerVal[a, b, c] =
  dFinite[thetaOf[a], thetaOf[b], thetaOf[c], ksqOf[a], ksqOf[b], ksqOf[c]];
BphiKerVal[a_, b_, c_] := BphiKerVal[a, b, c] =
  bFinite[thetaOf[a], thetaOf[b], thetaOf[c], ksqOf[a], ksqOf[b], ksqOf[c]];
H3kerVal[a_, b_, c_] :=
  1/4 (ksqOf[b] + ksqOf[c] - ksqOf[a] - 2 thetaOf[b] thetaOf[c]);

(* ------------------------------------------------------------------ *)
(*   Build Hamiltonians, W3, W4 for a given mode set                  *)
(* ------------------------------------------------------------------ *)

buildSystem[modeSet_] := Module[
  {H2loc, H3loc, W3loc, N2PhiLoc, S4H4loc, S4W3loc, S4Rawloc,
   toN, fromN, normalVarsLoc, H4loc,
   bracketLoc, splitFunc, solveFunc,
   W4H4loc, W4W3loc, W4Totloc},

  bracketLoc[F_, G_] := Sum[
    D[F, z[k]] D[G, p[-k]] - D[F, p[k]] D[G, z[-k]],
    {k, modeSet}];

  H2loc = 1/2 Sum[z[k] z[-k] + thetaOf[k] p[k] p[-k], {k, modeSet}];

  H3loc = Sum[
    If[MemberQ[modeSet, -(a+b)] && a+b != 0,
      z[-(a+b)] p[a] p[b] H3kerVal[-(a+b), a, b], 0],
    {a, modeSet}, {b, modeSet}];

  W3loc = Sum[
    If[a+b+c == 0,
      DkerVal[a,b,c] p[a] p[b] p[c] + BphiKerVal[a,b,c] z[a] z[b] p[c], 0],
    {a, modeSet}, {b, modeSet}, {c, modeSet}];

  N2PhiLoc[k_] := N2PhiLoc[k] = Sum[
    If[MemberQ[modeSet, d] && MemberQ[modeSet, k-b-d],
      z[b] z[k-b-d] p[d] (
        thetaOf[k] thetaOf[k-b] thetaOf[d]
        - 1/2 thetaOf[k] ksqOf[d]
        - 1/2 ksqOf[k] thetaOf[d]), 0],
    {b, modeSet}, {d, modeSet}];

  S4H4loc = 1/2 Sum[
    If[MemberQ[modeSet, -k], p[-k] N2PhiLoc[k], 0], {k, modeSet}];

  S4W3loc = Expand[
    -bracketLoc[H3loc, W3loc] + 1/2 bracketLoc[bracketLoc[H2loc, W3loc], W3loc]];

  S4Rawloc = Expand[S4H4loc + S4W3loc];

  toN = Flatten[Table[
    {z[k] -> (aa[k]+bb[k])/2,
     p[k] -> (aa[k]-bb[k])/(2 I omegaOf[k])},
    {k, modeSet}]];

  fromN = Flatten[Table[
    {aa[k] -> z[k] + I omegaOf[k] p[k],
     bb[k] -> z[k] - I omegaOf[k] p[k]},
    {k, modeSet}]];

  normalVarsLoc = Flatten[Table[{aa[k], bb[k]}, {k, modeSet}]];

  normalWeightLoc[mon_] := Module[{pa, pb},
    pa = Total[Table[Exponent[mon, aa[k]] omegaOf[k], {k, modeSet}]];
    pb = Total[Table[Exponent[mon, bb[k]] omegaOf[k], {k, modeSet}]];
    pa - pb];

  splitFunc[poly_] := Module[{terms, mon, wt, res, nonres},
    terms = If[Head[Expand[poly]] === Plus,
      List @@ Expand[poly], {Expand[poly]}];
    res = {}; nonres = {};
    Do[
      mon = Times @@ Table[var^Exponent[term, var], {var, normalVarsLoc}];
      wt = normalWeightLoc[mon];
      If[TrueQ[Simplify[wt == 0]],
        AppendTo[res, term], AppendTo[nonres, term]],
      {term, terms}];
    <|"nonres" -> Total[nonres], "res" -> Total[res]|>];

  solveFunc[poly_] := Module[{terms, mon, coeff, wt, solved},
    terms = If[Head[Expand[poly]] === Plus,
      List @@ Expand[poly], {Expand[poly]}];
    solved = Table[
      mon = Times @@ Table[var^Exponent[term, var], {var, normalVarsLoc}];
      coeff = Together[term/mon];
      wt = normalWeightLoc[mon];
      If[TrueQ[Simplify[wt == 0]], 0, -coeff mon/(I wt)],
      {term, terms}];
    Total[solved]];

  W4H4loc = Expand[
    solveFunc[splitFunc[Expand[S4H4loc /. toN]]["nonres"]] /. fromN];
  W4W3loc = Expand[
    solveFunc[splitFunc[Expand[S4W3loc /. toN]]["nonres"]] /. fromN];
  W4Totloc = Expand[W4H4loc + W4W3loc];

  <|"W3" -> W3loc, "W4Total" -> W4Totloc, "H2" -> H2loc,
    "bracketFn" -> bracketLoc|>
];

(* ------------------------------------------------------------------ *)
(*   Substitution rules                                                *)
(* ------------------------------------------------------------------ *)
rulesSigma6 = {
  T[1] -> sigma, T[2] -> 2 Tanh[2 ArcTanh[sigma]],
  T[3] -> 3 Tanh[3 ArcTanh[sigma]], T[4] -> 4 Tanh[4 ArcTanh[sigma]],
  T[5] -> 5 Tanh[5 ArcTanh[sigma]], T[6] -> 6 Tanh[6 ArcTanh[sigma]]};

sig6[expr_] := FullSimplify[FunctionExpand[expr /. rulesSigma6]];

stokesC33 = (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6);

(* extract cubic coeff of eps1^n1 eps2^n2 eps3^n3, scaled by 2 *)
extractC[expr_, rules_, n1_, n2_, n3_] :=
  2 Expand[Coefficient[Coefficient[Coefficient[
    expr /. rules, eps1, n1], eps2, n2], eps3, n3]];

biRules[modes_] := Flatten[Table[
  {z[k] -> Which[Abs[k]==1, eps1/2, Abs[k]==2, eps2/2, True, 0],
   p[k] -> Which[
     k == 1, -I eps1/(2 Sqrt[T[1]]), k == -1, I eps1/(2 Sqrt[T[1]]),
     k == 2, -I eps2/(2 Sqrt[T[2]]), k == -2, I eps2/(2 Sqrt[T[2]]),
     True, 0]},
  {k, modes}]];

triRules[modes_] := Flatten[Table[
  {z[k] -> Which[Abs[k]==1,eps1/2, Abs[k]==2,eps2/2, Abs[k]==3,eps3/2, True,0],
   p[k] -> Which[
     k==1, -I eps1/(2 Sqrt[T[1]]),  k==-1, I eps1/(2 Sqrt[T[1]]),
     k==2, -I eps2/(2 Sqrt[T[2]]),  k==-2, I eps2/(2 Sqrt[T[2]]),
     k==3, -I eps3/(2 Sqrt[T[3]]),  k==-3, I eps3/(2 Sqrt[T[3]]),
     True, 0]},
  {k, modes}]];

mapH3[m_, W3_] := Expand[D[W3, p[-m]] + 1/2 Expand[
  Sum[D[D[W3, p[-m]], z[k]] D[W3, p[-k]] -
      D[D[W3, p[-m]], p[k]] D[W3, z[-k]],
      {k, Range[1,6]}]]];

mapFull[m_, W3_, W4T_] := Expand[mapH3[m, W3] + D[W4T, p[-m]]];

(* ================================================================== *)
(*   PART A: mode set {-3,...,3}, cross terms at output mode 3        *)
(* ================================================================== *)
Print["Building system on modes {-3,...,3} ..."];
modesA = Flatten[{Range[-3,-1], Range[1,3]}];
sysA = buildSystem[modesA];
W3A  = sysA["W3"];
W4TA = sysA["W4Total"];
biA  = biRules[modesA];
triA = triRules[modesA];

Print[""];
Print["=== PART A: mode set {-3,..,3}, output mode 3 ==="];
Print[""];

Print["A1. Single-freq: (1,1,1) -> 3, eps1^3"];
mH3A = mapH3[3, W3A];  mFullA = mapFull[3, W3A, W4TA];
cA1h3   = sig6[extractC[mH3A, biA, 3, 0, 0]];
cA1full = sig6[extractC[mFullA, biA, 3, 0, 0]];
Print["  H3-only = ", cA1h3];
Print["  H3+H4   = ", cA1full];
Print["  Stokes  = ", stokesC33];
Print["  H3-only residual = ", sig6[cA1h3 - (3(1+sigma^2)(54+39sigma^2+70sigma^4-28sigma^6-4sigma^8-3sigma^10))/(64 sigma^6 (9+14sigma^2+9sigma^4))]];
Print["  H3+H4  residual  = ", sig6[cA1full - stokesC33]];
Print[""];

Print["A2. Bichromatic cross: (-1,2,2) -> 3, eps1*eps2^2"];
cA2h3   = sig6[extractC[mH3A, biA, 1, 2, 0]];
cA2full = sig6[extractC[mFullA, biA, 1, 2, 0]];
Print["  H3-only = ", cA2h3];
Print["  H3+H4   = ", cA2full];
Print["  Delta   = ", sig6[cA2full - cA2h3]];
Print["  DW DW limit sigma->1: H3-only = ", FullSimplify[cA2h3 /. sigma->1]];
Print["  DW DW limit sigma->1: H3+H4   = ", FullSimplify[cA2full /. sigma->1]];
Print[""];

Print["A3. Trichromatic cross: (-1,1,3) -> 3, eps1^2*eps3"];
cA3h3   = sig6[extractC[mH3A, triA, 2, 0, 1]];
cA3full = sig6[extractC[mapFull[3,W3A,W4TA], triA, 2, 0, 1]];
Print["  H3-only = ", cA3h3];
Print["  H3+H4   = ", cA3full];
Print["  Delta   = ", sig6[cA3full - cA3h3]];
Print["  DW limit sigma->1: H3-only = ", FullSimplify[cA3h3 /. sigma->1]];
Print["  DW limit sigma->1: H3+H4   = ", FullSimplify[cA3full /. sigma->1]];
Print[""];

Print["A4. Trichromatic cross: (-2,2,3) -> 3, eps2^2*eps3"];
cA4h3   = sig6[extractC[mH3A, triA, 0, 2, 1]];
cA4full = sig6[extractC[mapFull[3,W3A,W4TA], triA, 0, 2, 1]];
Print["  H3-only = ", cA4h3];
Print["  H3+H4   = ", cA4full];
Print["  Delta   = ", sig6[cA4full - cA4h3]];
Print["  DW limit sigma->1: H3-only = ", FullSimplify[cA4h3 /. sigma->1]];
Print["  DW limit sigma->1: H3+H4   = ", FullSimplify[cA4full /. sigma->1]];
Print[""];

(* ================================================================== *)
(*   PART B: mode set {1,...,6}, cross terms at output modes 4 5 6    *)
(* ================================================================== *)
Print["Building system on modes {1,...,6} ..."];
modesB = Range[1, 6];
sysB = buildSystem[modesB];
W3B  = sysB["W3"];
W4TB = sysB["W4Total"];
biB  = biRules[modesB];
triB = triRules[modesB];

rulesSigmaB = {
  T[1] -> sigma, T[2] -> 2 Tanh[2 ArcTanh[sigma]],
  T[3] -> 3 Tanh[3 ArcTanh[sigma]], T[4] -> 4 Tanh[4 ArcTanh[sigma]],
  T[5] -> 5 Tanh[5 ArcTanh[sigma]], T[6] -> 6 Tanh[6 ArcTanh[sigma]]};
sigB[expr_] := FullSimplify[FunctionExpand[expr /. rulesSigmaB]];

Print[""];
Print["=== PART B: mode set {-6,..,6}, output modes 4 5 6 ==="];
Print[""];

Print["B1. Bichromatic: (1,1,2) -> 4, eps1^2*eps2"];
mH3B4  = mapH3[4, W3B];   mFullB4 = mapFull[4, W3B, W4TB];
cB1h3   = sigB[extractC[mH3B4,  biB, 2, 1, 0]];
cB1full = sigB[extractC[mFullB4, biB, 2, 1, 0]];
Print["  H3-only = ", cB1h3];
Print["  H3+H4   = ", cB1full];
Print["  Delta   = ", sigB[cB1full - cB1h3]];
Print["  DW limit sigma->1: H3-only = ", FullSimplify[cB1h3 /. sigma->1]];
Print["  DW limit sigma->1: H3+H4   = ", FullSimplify[cB1full /. sigma->1]];
Print[""];

Print["B2. Bichromatic: (1,2,2) -> 5, eps1*eps2^2"];
mH3B5  = mapH3[5, W3B];   mFullB5 = mapFull[5, W3B, W4TB];
cB2h3   = sigB[extractC[mH3B5,  biB, 1, 2, 0]];
cB2full = sigB[extractC[mFullB5, biB, 1, 2, 0]];
Print["  H3-only = ", cB2h3];
Print["  H3+H4   = ", cB2full];
Print["  Delta   = ", sigB[cB2full - cB2h3]];
Print["  DW limit sigma->1: H3-only = ", FullSimplify[cB2h3 /. sigma->1]];
Print["  DW limit sigma->1: H3+H4   = ", FullSimplify[cB2full /. sigma->1]];
Print[""];

Print["B3. Single-freq: (2,2,2) -> 6, eps2^3"];
mH3B6  = mapH3[6, W3B];   mFullB6 = mapFull[6, W3B, W4TB];
cB3h3   = sigB[extractC[mH3B6,  biB, 0, 3, 0]];
cB3full = sigB[extractC[mFullB6, biB, 0, 3, 0]];
Print["  H3-only = ", cB3h3];
Print["  H3+H4   = ", cB3full];
Print["  DW limit sigma->1: H3-only = ", FullSimplify[cB3h3 /. sigma->1]];
Print["  DW limit sigma->1: H3+H4   = ", FullSimplify[cB3full /. sigma->1]];
Print[""];

Print["B4. Trichromatic: (1,1,3) -> 5, eps1^2*eps3"];
cB4h3   = sigB[extractC[mH3B5,  triB, 2, 0, 1]];
cB4full = sigB[extractC[mFullB5, triB, 2, 0, 1]];
Print["  H3-only = ", cB4h3];
Print["  H3+H4   = ", cB4full];
Print["  Delta   = ", sigB[cB4full - cB4h3]];
Print["  DW limit sigma->1: H3-only = ", FullSimplify[cB4h3 /. sigma->1]];
Print["  DW limit sigma->1: H3+H4   = ", FullSimplify[cB4full /. sigma->1]];
Print[""];

Print["B5. Trichromatic: (1,2,3) -> 6, eps1*eps2*eps3"];
cB5h3   = sigB[extractC[mH3B6,  triB, 1, 1, 1]];
cB5full = sigB[extractC[mFullB6, triB, 1, 1, 1]];
Print["  H3-only = ", cB5h3];
Print["  H3+H4   = ", cB5full];
Print["  Delta   = ", sigB[cB5full - cB5h3]];
Print["  DW limit sigma->1: H3-only = ", FullSimplify[cB5h3 /. sigma->1]];
Print["  DW limit sigma->1: H3+H4   = ", FullSimplify[cB5full /. sigma->1]];
Print[""];

Print["=== Summary ==="];
Print["Cross terms where H3+H4 Delta != 0 show where quartic absorption adds"];
Print["new observable content beyond the H3-only map."];
Print["Deep-water limit (sigma->1) should match known Stokes/VWA values."];
