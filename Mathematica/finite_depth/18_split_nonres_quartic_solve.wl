(* ::Package:: *)
(*  Solve the quartic homological equation only on non-resonant source pieces.

    This script separates the raw quartic source into:

      S4H4
      S4W3Induced = -{H3,W3} + 1/2 {{H2,W3},W3}
      S4Raw = S4H4 + S4W3Induced

    and solves

      {H2, W4H4}    = -S4H4Nonres
      {H2, W4W3}    = -S4W3Nonres
      {H2, W4Total} = -S4RawNonres

    The split-and-add step is only for the quartic source solve, not for the
    full canonical transformation.
*)

ClearAll["Global`*"];

$Assumptions = (And @@ Table[T[n] > 0, {n, 1, 6}]) &&
  Element[sigma, Reals] && 0 < sigma < 1;

modes = {-3, -2, -1, 1, 2, 3};

theta[0] := 0;
theta[n_Integer] /; n != 0 := T[Abs[n]];
omega[n_Integer] /; n != 0 := Sqrt[theta[n]];
ksq[n_Integer] := n^2;

denBase[t1_, t2_, t3_] :=
  t1 (t2 + t3 - t1) + t2 (t1 + t3 - t2) + t3 (t1 + t2 - t3);

dFinite[t1_, t2_, t3_, q1_, q2_, q3_] := Module[{den, num},
  den = denBase[t1, t2, t3];
  If[t1 === 0 || t2 === 0 || t3 === 0 || TrueQ[den == 0],
    0,
    num = q1 (t1^2 - (t2 - t3)^2)
        + q2 (t2^2 - (t1 - t3)^2)
        + q3 (t3^2 - (t1 - t2)^2)
        - 2 t1 t2 t3 (t1 + t2 + t3);
    num/(12 den)
  ]
];

bFinite[t1_, t2_, t3_, q1_, q2_, q3_] := Module[{den, num},
  den = denBase[t1, t2, t3];
  If[TrueQ[den == 0],
    0,
    num = t3 (t3 (t1 + t2) - t1^2 - t2^2)
        + t3 (q1 + q2 - 2 q3)
        + (t1 - t2) (q1 - q2);
    num/(2 den)
  ]
];

Dker[a_, b_, c_] := Dker[a, b, c] =
  dFinite[theta[a], theta[b], theta[c], ksq[a], ksq[b], ksq[c]];

BphiKer[a_, b_, c_] := BphiKer[a, b, c] =
  bFinite[theta[a], theta[b], theta[c], ksq[a], ksq[b], ksq[c]];

H3ker[a_, b_, c_] :=
  1/4 (ksq[b] + ksq[c] - ksq[a] - 2 theta[b] theta[c]);

bracketRaw[F_, G_] :=
  Sum[
    D[F, z[k]] D[G, p[-k]] - D[F, p[k]] D[G, z[-k]],
    {k, modes}
  ];

H2 = 1/2 Sum[z[k] z[-k] + theta[k] p[k] p[-k], {k, modes}];

H3 = Sum[
    If[MemberQ[modes, -(a + b)] && a + b != 0,
      z[-(a + b)] p[a] p[b] H3ker[-(a + b), a, b],
      0
    ],
    {a, modes}, {b, modes}
  ];

W3 = Sum[
    If[a + b + c == 0,
      Dker[a, b, c] p[a] p[b] p[c]
      + BphiKer[a, b, c] z[a] z[b] p[c],
      0
    ],
    {a, modes}, {b, modes}, {c, modes}
  ];

N2Phi[k_] := N2Phi[k] = Sum[
    If[MemberQ[modes, d] && MemberQ[modes, k - b - d],
      z[b] z[k - b - d] p[d] (
        theta[k] theta[k - b] theta[d]
        - 1/2 theta[k] ksq[d]
        - 1/2 ksq[k] theta[d]
      ),
      0
    ],
    {b, modes}, {d, modes}
  ];

S4H4 = 1/2 Sum[
    If[MemberQ[modes, -k], p[-k] N2Phi[k], 0],
    {k, modes}
  ];

S4W3Induced1Raw = Expand[bracketRaw[H3, W3]];
S4W3Induced2Raw = Expand[1/2 bracketRaw[bracketRaw[H2, W3], W3]];
S4W3Induced = Expand[-S4W3Induced1Raw + S4W3Induced2Raw];
S4Raw = Expand[S4H4 + S4W3Induced];

toNormalRules = Flatten[Table[
    {
      z[k] -> (aa[k] + bb[k])/2,
      p[k] -> (aa[k] - bb[k])/(2 I omega[k])
    },
    {k, modes}
  ]];

fromNormalRules = Flatten[Table[
    {
      aa[k] -> z[k] + I omega[k] p[k],
      bb[k] -> z[k] - I omega[k] p[k]
    },
    {k, modes}
  ]];

normalVars = Flatten[Table[{aa[k], bb[k]}, {k, modes}]];

normalWeight[mon_] := Module[{pa, pb},
  pa = Total[Table[Exponent[mon, aa[k]] omega[k], {k, modes}]];
  pb = Total[Table[Exponent[mon, bb[k]] omega[k], {k, modes}]];
  pa - pb
];

splitResonance[poly_] := Module[{terms, mon, weight, resonant, nonres},
  terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
  resonant = {};
  nonres = {};
  Do[
    mon = Times @@ Table[var^Exponent[term, var], {var, normalVars}];
    weight = normalWeight[mon];
    If[TrueQ[Simplify[weight == 0]],
      resonant = Append[resonant, term],
      nonres = Append[nonres, term]
    ],
    {term, terms}
  ];
  <|
    "resonant" -> Total[resonant],
    "nonres" -> Total[nonres],
    "resonantCount" -> Length[resonant],
    "nonresCount" -> Length[nonres]
  |>
];

solveHomologicalForNegativeSource[poly_] := Module[{terms, mon, coeff, weight, solved},
  terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
  solved = Table[
    mon = Times @@ Table[var^Exponent[term, var], {var, normalVars}];
    coeff = Together[term/mon];
    weight = normalWeight[mon];
    If[TrueQ[Simplify[weight == 0]], 0, -coeff mon/(I weight)],
    {term, terms}
  ];
  Total[solved]
];

monomialFromPowers[powers_] := Times @@ MapThread[#1^#2 &, {normalVars, powers}];

coefficientCheck[expr_] := Module[{rules, numericFiltered, checked, sigmaSamples},
  sigmaSamples = {sigma -> 1/5, sigma -> 2/3, sigma -> 9/10};
  rules = CoefficientRules[Expand[expr], normalVars];
  numericFiltered = Select[
    rules,
    Module[{vals},
      vals = Quiet[N[Table[#[[2]] /. rulesSigma /. s, {s, sigmaSamples}], 40]];
      Max[Abs[vals]] > 10^-20
    ] &
  ];
  checked = Table[
    <|
      "powers" -> rule[[1]],
      "monomial" -> monomialFromPowers[rule[[1]]],
      "coeff" -> Simplify[rule[[2]] /. rulesSigma]
    |>,
    {rule, numericFiltered}
  ];
  Select[checked, Not[TrueQ[#["coeff"] === 0]] &]
];

rulesSigma = {
  T[1] -> sigma,
  T[2] -> 2 Tanh[2 ArcTanh[sigma]],
  T[3] -> 3 Tanh[3 ArcTanh[sigma]],
  T[4] -> 4 Tanh[4 ArcTanh[sigma]],
  T[5] -> 5 Tanh[5 ArcTanh[sigma]],
  T[6] -> 6 Tanh[6 ArcTanh[sigma]]
};

numericStateRules = Flatten@Table[
   {
     aa[k] -> N[(1 + 0.3 k)/(10 + k^2), 30],
     bb[k] -> N[(1 - 0.2 k)/(12 + k^2), 30]
   },
   {k, modes}
];

sigmaSamples = {sigma -> 1/5, sigma -> 2/3, sigma -> 9/10};

polySummary[expr_, n_: 4] := Module[{terms},
  terms = If[TrueQ[Expand[expr] === 0], {}, List @@ Expand[expr]];
  <|
    "count" -> Length[terms],
    "sample" -> Take[terms, UpTo[n]]
  |>
];

S4H4Split = splitResonance[Expand[S4H4 /. toNormalRules]];
S4W3Split = splitResonance[Expand[S4W3Induced /. toNormalRules]];
S4RawSplit = splitResonance[Expand[S4Raw /. toNormalRules]];

S4H4Nonres = Expand[S4H4Split["nonres"]];
S4W3Nonres = Expand[S4W3Split["nonres"]];
S4RawNonres = Expand[S4RawSplit["nonres"]];

W4H4 = Expand[solveHomologicalForNegativeSource[S4H4Nonres] /. fromNormalRules];
W4W3 = Expand[solveHomologicalForNegativeSource[S4W3Nonres] /. fromNormalRules];
W4Total = Expand[W4H4 + W4W3];

residualH4 = Expand[(bracketRaw[H2, W4H4] + S4H4) /. toNormalRules];
residualW3 = Expand[(bracketRaw[H2, W4W3] + S4W3Induced) /. toNormalRules];
residualTotal = Expand[(bracketRaw[H2, W4Total] + S4Raw) /. toNormalRules];

residualH4Nonres = Expand[(bracketRaw[H2, W4H4] + S4H4Nonres /. fromNormalRules) /. toNormalRules];
residualW3Nonres = Expand[(bracketRaw[H2, W4W3] + S4W3Nonres /. fromNormalRules) /. toNormalRules];
residualTotalNonres = Expand[(bracketRaw[H2, W4Total] + S4RawNonres /. fromNormalRules) /. toNormalRules];

residualH4Samples = Table[Quiet[N[residualH4Nonres /. rulesSigma /. s /. numericStateRules, 30]], {s, sigmaSamples}];
residualW3Samples = Table[Quiet[N[residualW3Nonres /. rulesSigma /. s /. numericStateRules, 30]], {s, sigmaSamples}];
residualTotalSamples = Table[Quiet[N[residualTotalNonres /. rulesSigma /. s /. numericStateRules, 30]], {s, sigmaSamples}];

residualH4Terms = coefficientCheck[residualH4Nonres];
residualW3Terms = coefficientCheck[residualW3Nonres];
residualTotalTerms = coefficientCheck[residualTotalNonres];

Print["Split non-resonant quartic solve from frozen W3"];
Print[""];
Print["Non-resonant source counts {H4, W3-induced, raw}: ",
  {
    S4H4Split["nonresCount"],
    S4W3Split["nonresCount"],
    S4RawSplit["nonresCount"]
  }
];
Print["Resonant source counts {H4, W3-induced, raw}: ",
  {
    S4H4Split["resonantCount"],
    S4W3Split["resonantCount"],
    S4RawSplit["resonantCount"]
  }
];
Print[""];
Print["Named non-resonant source pieces:"];
Print["  S4H4Nonres term count = ", polySummary[S4H4Nonres]["count"]];
Print[""];
Print["  S4W3Nonres term count = ", polySummary[S4W3Nonres]["count"]];
Print[""];
Print["  S4RawNonres term count = ", polySummary[S4RawNonres]["count"]];
Print[""];
Print["Named quartic generators constructed:"];
Print["  W4H4 ready for {H2,W4H4} = -S4H4Nonres"];
Print["  W4W3 ready for {H2,W4W3} = -S4W3Nonres"];
Print["  W4Total = W4H4 + W4W3 ready for {H2,W4Total} = -S4RawNonres"];
Print[""];
Print["Sampled residual checks for {H2,W4} + source = 0:"];
Print["  W4H4 samples   = ", residualH4Samples];
Print["  W4W3 samples   = ", residualW3Samples];
Print["  W4Total samples = ", residualTotalSamples];
Print[""];
Print["Exact coefficient residual counts:"];
Print["  W4H4 residual term count   = ", Length[residualH4Terms]];
Print["  W4W3 residual term count   = ", Length[residualW3Terms]];
Print["  W4Total residual term count = ", Length[residualTotalTerms]];
Print[""];
Print["Residual samples (up to 6 terms each) if nonzero:"];
Print["  W4H4   -> ", Take[residualH4Terms, UpTo[6]]];
Print["  W4W3   -> ", Take[residualW3Terms, UpTo[6]]];
Print["  W4Total -> ", Take[residualTotalTerms, UpTo[6]]];
