(* ::Package:: *)
(*
  Direct H4 convolution and overall H3+H4 absorption check.

  This script addresses two implementation risks from the previous attempts:

    1. The quartic Hamiltonian H4 should be built directly from the operator
       formula

         N2 psi = G0 eta G0 eta G0 psi
                  + 1/2 G0(eta^2 Delta psi)
                  + 1/2 Delta(eta^2 G0 psi)

       rather than by trusting a hand-written ordered four-wave kernel.

    2. The Stokes comparison should use the overall coordinate contribution

         1/2 {{zeta,W3},W3} + {zeta,W4}

       while W4 is obtained only from the non-resonant part of K4.

  The "finite-depth gauge" check below is applied at the observable level:
  delta_W4(sigma) -> delta_W4(sigma) - delta_W4(1), so the quartic correction
  cannot change the known deep-water C33=3/8 result.
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

Dker[dest_Integer, p1_Integer, r_Integer] :=
  Dker[dest, p1, r] =
    dFinite[theta[dest], theta[p1], theta[r], ksq[dest], ksq[p1], ksq[r]];

BphiKer[dest_Integer, p1_Integer, r_Integer] :=
  BphiKer[dest, p1, r] =
    bFinite[theta[dest], theta[p1], theta[r], ksq[dest], ksq[p1], ksq[r]];

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

(* Direct convolution from N2.  Mode k is the output mode of N2[phi]. *)
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

H4OrderedKer[a_, b_, c_, d_] :=
  1/2 theta[a] theta[c + d] theta[d]
  - 1/4 (theta[a] ksq[d] + ksq[a] theta[d]);

H4Ordered = Sum[
    If[a + b + c + d == 0,
      z[b] z[c] p[a] p[d] H4OrderedKer[a, b, c, d],
      0
    ],
    {a, modes}, {b, modes}, {c, modes}, {d, modes}
  ];

h4Difference = Together[Expand[H4Direct - H4Ordered]];

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
    If[TrueQ[FullSimplify[weight == 0]],
      resonant = Append[resonant, term],
      nonres = Append[nonres, term]
    ],
    {term, terms}
  ];
  <|"resonant" -> Total[resonant], "nonres" -> Total[nonres],
    "resonantCount" -> Length[resonant], "nonresCount" -> Length[nonres]|>
];

solveHomologicalNonres[poly_] := Module[{terms, mon, coeff, weight, solved},
  terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
  solved = Table[
    mon = Times @@ Table[var^Exponent[term, var], {var, normalVars}];
    coeff = Together[term/mon];
    weight = normalWeight[mon];
    If[TrueQ[FullSimplify[weight == 0]], 0, coeff mon/(I weight)],
    {term, terms}
  ];
  Total[solved]
];

K4Direct = Expand[H4Direct + 1/2 bracketRaw[H3, W3]];
K4Normal = Expand[K4Direct /. toNormalRules];
K4Split = splitResonance[K4Normal];
W4Normal = solveHomologicalNonres[K4Normal];
W4Back = Expand[W4Normal /. fromNormalRules];

computeDeltaW4[bracketSign_, inverseSign_] := Module[
  {K4Local, W4NormalLocal, W4BackLocal},
  K4Local = Expand[(H4Direct + bracketSign 1/2 bracketRaw[H3, W3]) /. toNormalRules];
  W4NormalLocal = inverseSign solveHomologicalNonres[K4Local];
  W4BackLocal = Expand[W4NormalLocal /. fromNormalRules];
  cosC33FromPositiveMode[bracketRaw[z[3], W4BackLocal]]
];

linearWaveRules = Flatten[Table[
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

rulesSigma = {
  T[1] -> sigma,
  T[2] -> 2 Tanh[2 ArcTanh[sigma]],
  T[3] -> 3 Tanh[3 ArcTanh[sigma]],
  T[4] -> 4 Tanh[4 ArcTanh[sigma]],
  T[5] -> 5 Tanh[5 ArcTanh[sigma]],
  T[6] -> 6 Tanh[6 ArcTanh[sigma]]
};

cosC33FromPositiveMode[expr_] :=
  FullSimplify[2 Coefficient[Expand[expr /. linearWaveRules], eps, 3] /. rulesSigma];

h3OnlyC33 = cosC33FromPositiveMode[
  1/2 bracketRaw[bracketRaw[z[3], W3], W3]
];

deltaW4C33 = cosC33FromPositiveMode[bracketRaw[z[3], W4Back]];
deltaW4C33Deep = FullSimplify[deltaW4C33 /. sigma -> 1];
deltaW4C33FiniteDepthGauge = FullSimplify[deltaW4C33 - deltaW4C33Deep];

stokesC33 = (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6);
requiredDeltaC33 = FullSimplify[stokesC33 - h3OnlyC33];

overallC33 = FullSimplify[h3OnlyC33 + deltaW4C33];
overallGaugeC33 = FullSimplify[h3OnlyC33 + deltaW4C33FiniteDepthGauge];

variantData = Flatten[
  Table[
    delta = computeDeltaW4[bs, is];
    deltaGauge = FullSimplify[delta - (delta /. sigma -> 1)];
    <|
      "bracketSign" -> bs,
      "inverseSign" -> is,
      "s4" -> s4,
      "delta" -> FullSimplify[s4 delta],
      "deltaGauge" -> FullSimplify[s4 deltaGauge],
      "residual" -> FullSimplify[h3OnlyC33 + s4 delta - stokesC33],
      "residualGauge" -> FullSimplify[h3OnlyC33 + s4 deltaGauge - stokesC33],
      "deepGauge" -> FullSimplify[(h3OnlyC33 + s4 deltaGauge) /. sigma -> 1]
    |>,
    {bs, {-1, 1}}, {is, {-1, 1}}, {s4, {-1, 1}}
  ],
  2
];

zeroResidualVariants = Select[variantData, FullSimplify[#["residual"] == 0] &];
zeroGaugeResidualVariants = Select[variantData, FullSimplify[#["residualGauge"] == 0] &];

Print["Direct H4 convolution overall check"];
Print[""];
Print["H4 direct convolution equals previous ordered kernel? ",
  TrueQ[h4Difference === 0]
];
Print["  Number of terms in H4Direct-H4Ordered after expansion: ",
  If[h4Difference === 0, 0, Length[List @@ Expand[h4Difference]]]
];
Print[""];
Print["K4 resonance split in normal variables:"];
Print["  resonant term count:     ", K4Split["resonantCount"]];
Print["  non-resonant term count: ", K4Split["nonresCount"]];
Print[""];
Print["Single-frequency C33 comparison:"];
Print["  H3-only C33:                  ", h3OnlyC33];
Print["  Stokes/VWA C33:               ", stokesC33];
Print["  required H4-level delta:      ", requiredDeltaC33];
Print["  nonres W4 delta:              ", deltaW4C33];
Print["  nonres W4 delta deep-water:   ", deltaW4C33Deep];
Print["  finite-depth-gauge W4 delta:  ", deltaW4C33FiniteDepthGauge];
Print[""];
Print["Residuals:"];
Print["  H3-only - Stokes:             ", FullSimplify[h3OnlyC33 - stokesC33]];
Print["  H3+W4 - Stokes:               ", FullSimplify[overallC33 - stokesC33]];
Print["  H3+(W4-W4deep) - Stokes:      ", FullSimplify[overallGaugeC33 - stokesC33]];
Print[""];
Print["Focused sign scan using direct H4:"];
Do[
  Print[
    "  {bracketSign,inverseSign,s4} = ",
    {row["bracketSign"], row["inverseSign"], row["s4"]},
    " -> residual = ", row["residual"],
    "; gauge residual = ", row["residualGauge"],
    "; gauge deep-water C33 = ", row["deepGauge"]
  ],
  {row, variantData}
];
Print["  exact zero residual variants: ", zeroResidualVariants];
Print["  exact zero gauge residual variants: ", zeroGaugeResidualVariants];
Print[""];
Print["Deep-water checks sigma -> 1:"];
Print["  {H3-only, W4 delta, gauge delta, H3+gauge, Stokes}: ",
  FullSimplify[
    {h3OnlyC33, deltaW4C33, deltaW4C33FiniteDepthGauge,
      overallGaugeC33, stokesC33} /. sigma -> 1
  ]
];
