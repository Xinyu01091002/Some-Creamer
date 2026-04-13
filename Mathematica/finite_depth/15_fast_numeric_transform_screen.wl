(* ::Package:: *)
(*  Fast numeric screen for finite-depth H3/H4 transform conventions.

    This is a lighter companion to 13/14.  It keeps the exact
    single-frequency checks, but replaces the full quartic symbolic residual
    test by sampled Hamiltonian-residual evaluations in normal variables.
    The purpose is to iterate quickly and identify which transform
    convention is worth the heavier exact check.
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

H4Direct = 1/2 Sum[
    If[MemberQ[modes, -k], p[-k] N2Phi[k], 0],
    {k, modes}
  ];

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

solveHomologicalNonres[poly_] := Module[{terms, mon, coeff, weight, solved},
  terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
  solved = Table[
    mon = Times @@ Table[var^Exponent[term, var], {var, normalVars}];
    coeff = Together[term/mon];
    weight = normalWeight[mon];
    If[TrueQ[Simplify[weight == 0]], 0, coeff mon/(I weight)],
    {term, terms}
  ];
  Total[solved]
];

rulesSigma = {
  T[1] -> sigma,
  T[2] -> 2 Tanh[2 ArcTanh[sigma]],
  T[3] -> 3 Tanh[3 ArcTanh[sigma]],
  T[4] -> 4 Tanh[4 ArcTanh[sigma]],
  T[5] -> 5 Tanh[5 ArcTanh[sigma]],
  T[6] -> 6 Tanh[6 ArcTanh[sigma]]
};

singleFrequencyRules = Flatten[Table[
    {
      z[k] -> Which[k == 1, eps/2, k == -1, eps/2, True, 0],
      p[k] -> Which[
        k == 1, -I eps/(2 Sqrt[T[1]]),
        k == -1, I eps/(2 Sqrt[T[1]]),
        True, 0
      ]
    },
    {k, modes}
  ]];

exprToCubic[expr_, rules_] := Expand[Coefficient[Expand[expr /. rules], eps, 3]];

stokesC33 = (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6);
h3OnlyTarget =
  (3 (1 + sigma^2) (54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 - 4 sigma^8 - 3 sigma^10))/
    (64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

z3W3First = Expand[D[W3, p[-3]]];
z3W3Second = Expand[bracketRaw[z3W3First, W3]];
z3H3Map[s3_] := z[3] + s3 z3W3First + 1/2 s3^2 z3W3Second;

s3Data = Table[
   <|
     "s3" -> s3,
     "c33" -> Simplify[2 Coefficient[exprToCubic[z3H3Map[s3], singleFrequencyRules] /. rulesSigma, eps, 0]],
     "residual" -> Simplify[
       2 Coefficient[exprToCubic[z3H3Map[s3], singleFrequencyRules] /. rulesSigma, eps, 0] - h3OnlyTarget]
   |>,
   {s3, {-1, 1}}
];
chosenS3 = If[Length[Select[s3Data, TrueQ[#["residual"] === 0] &]] > 0,
  First[Select[s3Data, TrueQ[#["residual"] === 0] &]]["s3"],
  1
];

numericStateRules = Flatten@Table[
   {
     aa[k] -> N[(1 + 0.3 k)/(10 + k^2), 30],
     bb[k] -> N[(1 - 0.2 k)/(12 + k^2), 30]
   },
   {k, modes}
];

sigmaSamples = {sigma -> 1/5, sigma -> 2/3, sigma -> 9/10};

candidateSeeds = {
  <|"bracketSign" -> -1, "inverseSign" -> -1, "s4" -> 1|>,
  <|"bracketSign" -> -1, "inverseSign" -> 1, "s4" -> -1|>,
  <|"bracketSign" -> 1, "inverseSign" -> -1, "s4" -> 1|>,
  <|"bracketSign" -> 1, "inverseSign" -> 1, "s4" -> -1|>
};

candidateData = Table[
   Module[{seed, k4Normal, split, w4Normal, w4Back, adH2W4, r4Expr, r4Samples,
     z3Map, c33, c33Residual},
     seed = candidateSeeds[[idx]];
     k4Normal = Expand[(H4Direct + seed["bracketSign"] 1/2 bracketRaw[H3, W3]) /. toNormalRules];
     split = splitResonance[k4Normal];
     w4Normal = seed["inverseSign"] solveHomologicalNonres[split["nonres"]];
     w4Back = Expand[w4Normal /. fromNormalRules];
     adH2W4 = Expand[(bracketRaw[H2, w4Back]) /. toNormalRules];
     r4Expr = Expand[split["nonres"] - adH2W4];
     r4Samples = Table[
       N[r4Expr /. rulesSigma /. s /. numericStateRules, 30],
       {s, sigmaSamples}
     ];
     z3Map = z3H3Map[chosenS3] + seed["s4"] Expand[D[w4Back, p[-3]]];
     c33 = Simplify[2 Coefficient[exprToCubic[z3Map, singleFrequencyRules] /. rulesSigma, eps, 0]];
     c33Residual = Simplify[c33 - stokesC33];
     <|
       "bracketSign" -> seed["bracketSign"],
       "inverseSign" -> seed["inverseSign"],
       "s4" -> seed["s4"],
       "K4ResCount" -> split["resonantCount"],
       "K4NonresCount" -> split["nonresCount"],
       "singleResidual" -> c33Residual,
       "r4NumericSamples" -> r4Samples,
       "r4NumericMaxAbs" -> Max[Abs[r4Samples]]
     |>
   ],
   {idx, Length[candidateSeeds]}
];

Print["Fast numeric finite-depth H3/H4 transform screen"];
Print[""];
Print["Chosen W3 observable sign s3 from H3-only single-frequency check: ", chosenS3];
Print[""];
Do[
  Print[
    "  candidate {bracketSign,inverseSign,s4} = ",
    {row["bracketSign"], row["inverseSign"], row["s4"]},
    " -> K4{res,nonres} = ", {row["K4ResCount"], row["K4NonresCount"]},
    "; single residual = ", row["singleResidual"],
    "; numeric R4 samples = ", row["r4NumericSamples"],
    "; max abs sample = ", row["r4NumericMaxAbs"]
  ],
  {row, candidateData}
];
