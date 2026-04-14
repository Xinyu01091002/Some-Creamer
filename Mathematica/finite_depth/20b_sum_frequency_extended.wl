(* ::Package:: *)
(*  Sum-frequency cross-term coefficients, extended mode set {+-1,...,+-6}.
    Computes H3-only and H3+H4 cubic coefficients for:
      Bichromatic  {1,2}: outputs at modes 3, 4, 5, 6
      Trichromatic {1,2,3}: additional outputs at modes 5, 6

    FOCUS: sum-frequency triads only (a>0, b>0, c>0, output = a+b+c).
    Difference-frequency triads are excluded from this script.

    Strategy
    --------
    - Use numeric T[n] = n*tanh(n*kh) from the start; avoids symbolic bloat.
    - Use K4 = H4 - 1/2 {H3,W3} directly (this equals S4Raw because
      {H2,W3} = H3 collapses the double bracket: -[H3,W3] + 1/2[[H2,W3],W3]
      = -[H3,W3] + 1/2[H3,W3] = -1/2[H3,W3]).
    - Avoid FullSimplify; evaluate each result numerically.
    - Rebuild the full system for each kh value.

    Output convention
    -----------------
    For a cosine wave of physical amplitude A at mode k, set A[k] = A.
    The coefficient printed is 2*C where eta_nl[m] ~ (C/2) * prod(A[k_i]) * cos(m*x).
    This matches the Stokes C33 convention used in the single-frequency scripts.
*)

ClearAll["Global`*"];

(* ------------------------------------------------------------------ *)
(*  Core numerical system builder                                       *)
(* ------------------------------------------------------------------ *)

buildSystem[kh_] := Module[
  {modes, tv, tFn, oFn, ksqFn,
   dFiniteFn, bFiniteFn, DkerFn, BphiKerFn, H3kerFn,
   H2, H3, W3, H4, bracketH3W3, K4,
   toN, fromN, normalVarsLoc, nwtFn, splitNonresFn, solveFn,
   W4H4, W4W3, W4Total},

  modes = Flatten[{Range[-6, -1], Range[1, 6]}];

  (* numeric theta(n) = n*tanh(n*kh) *)
  tFn[n_] := If[n == 0, 0,
    If[kh === Infinity, N[Abs[n], 20], N[Abs[n] * Tanh[Abs[n] * kh], 20]]];
  oFn[n_] := If[n == 0, 0, Sqrt[tFn[n]]];
  ksqFn[n_] := n^2;

  dFiniteFn[t1_, t2_, t3_, q1_, q2_, q3_] := Module[{den, num},
    den = t1 (t2 + t3 - t1) + t2 (t1 + t3 - t2) + t3 (t1 + t2 - t3);
    If[t1 == 0 || t2 == 0 || t3 == 0 || den == 0, 0,
      num = q1 (t1^2-(t2-t3)^2) + q2 (t2^2-(t1-t3)^2) + q3 (t3^2-(t1-t2)^2)
          - 2 t1 t2 t3 (t1+t2+t3);
      num/(12 den)]];

  bFiniteFn[t1_, t2_, t3_, q1_, q2_, q3_] := Module[{den, num},
    den = t1 (t2 + t3 - t1) + t2 (t1 + t3 - t2) + t3 (t1 + t2 - t3);
    If[den == 0, 0,
      num = t3 (t3 (t1+t2) - t1^2 - t2^2) + t3 (q1+q2-2 q3) + (t1-t2)(q1-q2);
      num/(2 den)]];

  DkerFn[a_, b_, c_] := DkerFn[a, b, c] =
    dFiniteFn[tFn[a], tFn[b], tFn[c], ksqFn[a], ksqFn[b], ksqFn[c]];
  BphiKerFn[a_, b_, c_] := BphiKerFn[a, b, c] =
    bFiniteFn[tFn[a], tFn[b], tFn[c], ksqFn[a], ksqFn[b], ksqFn[c]];
  H3kerFn[a_, b_, c_] :=
    1/4 (ksqFn[b] + ksqFn[c] - ksqFn[a] - 2 tFn[b] tFn[c]);

  (* Poisson bracket on this mode set *)
  bracketFn[F_, G_] := Sum[
    D[F, z[k]] D[G, p[-k]] - D[F, p[k]] D[G, z[-k]], {k, modes}];

  H2 = 1/2 Sum[z[k] z[-k] + tFn[k] p[k] p[-k], {k, modes}];

  H3 = Expand[Sum[
    If[MemberQ[modes, -(a+b)] && a+b != 0,
      z[-(a+b)] p[a] p[b] H3kerFn[-(a+b), a, b], 0],
    {a, modes}, {b, modes}]];

  W3 = Expand[Sum[
    If[a + b + c == 0,
      DkerFn[a,b,c] p[a] p[b] p[c] + BphiKerFn[a,b,c] z[a] z[b] p[c], 0],
    {a, modes}, {b, modes}, {c, modes}]];

  (* H4 from N2 convolution: N2Phi[k] = sum_{b,d in modes, k-b-d in modes}
     z[b] z[k-b-d] p[d] * (theta_k * theta_{k-b} * theta_d
       - 1/2 theta_k * kd^2 - 1/2 kk^2 * theta_d) *)
  H4 = Expand[1/2 Sum[
    If[MemberQ[modes, -k],
      p[-k] Sum[
        If[MemberQ[modes, d] && MemberQ[modes, k-b-d],
          z[b] z[k-b-d] p[d] (
            tFn[k] tFn[k-b] tFn[d]
            - 1/2 tFn[k] ksqFn[d]
            - 1/2 ksqFn[k] tFn[d]), 0],
        {b, modes}, {d, modes}], 0],
    {k, modes}]];

  (* K4 = H4 - 1/2 {H3, W3}  (uses {H2,W3}=H3 to collapse double bracket) *)
  bracketH3W3 = Expand[bracketFn[H3, W3]];
  K4 = Expand[H4 - 1/2 bracketH3W3];

  (* Normal variables: aa[k] = z[k] + i omega[k] p[k],
                       bb[k] = z[k] - i omega[k] p[k]   *)
  toN = Flatten[Table[{
    z[k] -> (aa[k] + bb[k])/2,
    p[k] -> (aa[k] - bb[k])/(2 I oFn[k])}, {k, modes}]];
  fromN = Flatten[Table[{
    aa[k] -> z[k] + I oFn[k] p[k],
    bb[k] -> z[k] - I oFn[k] p[k]}, {k, modes}]];
  normalVarsLoc = Flatten[Table[{aa[k], bb[k]}, {k, modes}]];

  nwtFn[mon_] := Module[{pa, pb},
    pa = Total[Table[Exponent[mon, aa[k]] oFn[k], {k, modes}]];
    pb = Total[Table[Exponent[mon, bb[k]] oFn[k], {k, modes}]];
    pa - pb];

  splitNonresFn[poly_] := Module[{terms, mon, wt, nonres = {}},
    terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
    Do[
      mon = Times @@ Table[var^Exponent[t, var], {var, normalVarsLoc}];
      wt = nwtFn[mon];
      If[!TrueQ[Simplify[wt == 0, Assumptions -> True]], AppendTo[nonres, t]],
      {t, terms}];
    Total[nonres]];

  solveFn[poly_] := Module[{terms, mon, c, wt, solved},
    terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
    solved = Table[
      mon = Times @@ Table[var^Exponent[t, var], {var, normalVarsLoc}];
      c = Together[t/mon];
      wt = nwtFn[mon];
      If[TrueQ[Simplify[wt == 0, Assumptions -> True]], 0, -c mon/(I wt)],
      {t, terms}];
    Total[solved]];

  K4Normal = Expand[K4 /. toN];
  K4Nonres = splitNonresFn[K4Normal];
  W4Total = Expand[solveFn[K4Nonres] /. fromN];

  <| "modes" -> modes, "W3" -> W3, "W4Total" -> W4Total,
     "tFn" -> tFn, "oFn" -> oFn |>
];

(* ------------------------------------------------------------------ *)
(*  Coordinate map at output mode m (linear term z[m] assumed zero)   *)
(* ------------------------------------------------------------------ *)

mapH3at[m_, sys_] := Module[{W3 = sys["W3"], modes = sys["modes"]},
  Expand[
    D[W3, p[-m]] +
    1/2 Sum[
      D[D[W3, p[-m]], z[k]] D[W3, p[-k]] - D[D[W3, p[-m]], p[k]] D[W3, z[-k]],
      {k, modes}]]];

mapFullat[m_, sys_] :=
  Expand[mapH3at[m, sys] + D[sys["W4Total"], p[-m]]];

(* ------------------------------------------------------------------ *)
(*  Wave state substitution and coefficient extraction                 *)
(* ------------------------------------------------------------------ *)

(* Cosine wave: z[|k|] = z[-|k|] = A[|k|]/2,
                p[|k|]  = -i A[|k|] / (2 sqrt(theta(|k|))),
                p[-|k|] = +i A[|k|] / (2 sqrt(theta(|k|))) *)
waveSub[sys_, activeInputModes_] := Flatten[Table[
  Module[{o = sys["oFn"][k]},
    If[MemberQ[activeInputModes, k],
      {z[k] -> A[k]/2, z[-k] -> A[k]/2,
       p[k] -> -I A[k]/(2 o), p[-k] -> I A[k]/(2 o)},
      {z[k] -> 0, z[-k] -> 0, p[k] -> 0, p[-k] -> 0}]],
  {k, 1, 6}]];

(* extract cubic coefficient of A[1]^n1 * A[2]^n2 * A[3]^n3, scaled by 2 *)
getC[expr_, n1_, n2_, n3_] :=
  2 Coefficient[Coefficient[Coefficient[expr, A[1], n1], A[2], n2], A[3], n3];

(* ------------------------------------------------------------------ *)
(*  Sum-frequency triads to evaluate (output mode, {n1,n2,n3})        *)
(*  Only a>0, b>0, c>0 interactions: all mode indices positive.       *)
(* ------------------------------------------------------------------ *)

(* Bichromatic {1,2}: output modes 3,4,5,6 with powers of A[1] and A[2] *)
biTriads = {
  (* m=3 *) {3, {3,0,0}}, (* (1,1,1)->3   eps1^3 single-freq sanity *)
  (* m=4 *) {4, {2,1,0}}, (* (1,1,2)->4   eps1^2 eps2 *)
  (* m=5 *) {5, {1,2,0}}, (* (1,2,2)->5   eps1 eps2^2 *)
  (* m=6 *) {6, {0,3,0}}  (* (2,2,2)->6   eps2^3 single-freq for mode 2 *)
};

(* Trichromatic {1,2,3}: additional cross terms *)
triTriads = {
  (* m=5 *) {5, {2,0,1}}, (* (1,1,3)->5   eps1^2 eps3 *)
  (* m=6 *) {6, {1,1,1}}, (* (1,2,3)->6   eps1 eps2 eps3 *)
  (* m=5 *) {5, {0,1,2}}, (* (2,3,?)  wait: 2+3=5 needs third=0  ->
                              Actually (2,3) is only bichromatic {2,3},
                              omit; only write triads that are achievable *)
  (* m=6 *) {6, {0,0,3}}  (* (3,3,3)->6? No: 3+3=6 needs two, plus one=0?
                              Actually (2,2,2)=6 is above; (3,3,0) not sum-freq.
                              Add (1,2,3)->6 already above. *)
  {5, {0, 1, 2}},  (* (2,3,?): skip, 2+3=5 is only bichromatic from {2,3} *)
  {6, {3, 0, 0}},  (* (1,1,1,1,1,1)? No, these are cubic. (1,1,4)->6 is bichromatic
                      from {1,4}. For trichromatic {1,2,3}: (3,3,0) not sum-freq.
                      Actually (1,2,3)->6 is already listed. *)
  {7, {1, 0, 2}}   (* (1,3,3)->7  eps1 eps3^2  -- outside mode 6 range, skip *)
};

(* Clean trichromatic triad list (only achievable sum-freq within modes 1-6) *)
triTriads = {
  {5, {2, 0, 1}},   (* (1,1,3) -> 5 *)
  {6, {1, 1, 1}},   (* (1,2,3) -> 6 *)
  {6, {0, 0, 3}}    (* (3,3,0)? No, (3,3) is output 6, but needs zero = mode 0,
                       excluded. Keep only (1,2,3)->6 *)
};
triTriads = {
  {5, {2, 0, 1}},   (* (1,1,3) -> 5, eps1^2 eps3 *)
  {6, {1, 1, 1}}    (* (1,2,3) -> 6, eps1 eps2 eps3 *)
};

(* ------------------------------------------------------------------ *)
(*  Main computation loop over kh values                               *)
(* ------------------------------------------------------------------ *)

khValues = {1/2, 1, 2, Infinity};
khLabels = {"kh=0.5", "kh=1", "kh=2", "DW (kh=inf)"};

(* Reference: Stokes C33 = (27-9s^2+9s^4-3s^6)/(64 s^6) at sigma=tanh(kh) *)
stokesC33num[kh_] := If[kh === Infinity, 3/8,
  Module[{s = Tanh[kh]},
    N[(27-9s^2+9s^4-3s^6)/(64 s^6), 12]]];

Print["=== Sum-frequency cross terms: H3-only vs H3+H4 ==="];
Print["Mode set: {+-1, ..., +-6}.  All T[n] numeric."];
Print["Coefficient 2C where eta[m] ~ (C/2) * A[1]^n1 * A[2]^n2 * A[3]^n3 * cos(m x)"];
Print[""];

Do[
  kh = khValues[[idx]];
  label = khLabels[[idx]];
  Print[">>> ", label, " <<<"];
  Print["  Building system (W3 + W4Total) on modes {+-1,...,+-6} ..."];
  sys = buildSystem[kh];

  Print["  Computing maps at output modes 3,4,5,6 ..."];

  (* pre-build maps *)
  mH3  = <||>;  mFull = <||>;
  Do[
    mH3[m]   = mapH3at[m, sys];
    mFull[m] = mapFullat[m, sys],
    {m, {3, 4, 5, 6}}];

  (* bichromatic wave state {1,2} *)
  biSub = waveSub[sys, {1, 2}];
  biH3   = <||>;  biFull = <||>;
  Do[
    biH3[m]   = Expand[mH3[m]   /. biSub];
    biFull[m] = Expand[mFull[m] /. biSub],
    {m, {3, 4, 5, 6}}];

  (* trichromatic wave state {1,2,3} *)
  triSub = waveSub[sys, {1, 2, 3}];
  triH3  = <||>;  triFull = <||>;
  Do[
    triH3[m]   = Expand[mH3[m]   /. triSub];
    triFull[m] = Expand[mFull[m] /. triSub],
    {m, {5, 6}}];

  Print["  --- Bichromatic {1,2} ---"];
  Do[
    {m, pows} = triad;
    {n1, n2, n3} = pows;
    ch3   = Re[N[getC[biH3[m],   n1, n2, n3], 10]];
    cfull = Re[N[getC[biFull[m], n1, n2, n3], 10]];
    delta = cfull - ch3;
    ratio = If[Abs[ch3] > 10^-12, N[cfull/ch3, 8], "---"];
    Print["  (" <> ToString[n1] <> "," <> ToString[n2] <> "," <> ToString[n3] <>
          ")->" <> ToString[m] <> " :"];
    Print["    H3-only = " <> ToString[NumberForm[ch3, {10, 6}]]];
    Print["    H3+H4   = " <> ToString[NumberForm[cfull, {10, 6}]]];
    Print["    Delta   = " <> ToString[NumberForm[delta, {8, 5}]] <>
          "  ratio = " <> ToString[NumberForm[ratio, {8, 6}]]],
    {triad, biTriads}];

  (* single-freq sanity check *)
  Print["  [Stokes C33 ref: " <> ToString[N[stokesC33num[kh], 8]] <> "]"];
  Print[""];

  Print["  --- Trichromatic {1,2,3} ---"];
  Do[
    {m, pows} = triad;
    {n1, n2, n3} = pows;
    ch3   = Re[N[getC[triH3[m],   n1, n2, n3], 10]];
    cfull = Re[N[getC[triFull[m], n1, n2, n3], 10]];
    delta = cfull - ch3;
    ratio = If[Abs[ch3] > 10^-12, N[cfull/ch3, 8], "---"];
    Print["  (" <> ToString[n1] <> "," <> ToString[n2] <> "," <> ToString[n3] <>
          ")->" <> ToString[m] <> " :"];
    Print["    H3-only = " <> ToString[NumberForm[ch3, {10, 6}]]];
    Print["    H3+H4   = " <> ToString[NumberForm[cfull, {10, 6}]]];
    Print["    Delta   = " <> ToString[NumberForm[delta, {8, 5}]] <>
          "  ratio = " <> ToString[NumberForm[ratio, {8, 6}]]],
    {triad, triTriads}];
  Print[""],
  {idx, 1, Length[khValues]}];

Print["=== Done ==="];
Print[""];
Print["Interpretation guide"];
Print["  ratio > 1  : H3+H4 amplifies the coefficient relative to H3-only"];
Print["  ratio < 1  : H3+H4 reduces it"];
Print["  ratio = 1  : no quartic correction at this kh (deep-water case expected)"];
Print["  Compare sign and magnitude against MF12 numerical output"];
Print["  to assess whether H3+H4 moves each sum-frequency coefficient"];
Print["  in the correct direction."];
