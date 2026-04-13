(* ::Package:: *)
(*  Hamiltonian-first quartic validation from the frozen W3 basis.

    This is the formal mainline after 16-18:

      1. verify the frozen cubic W3 basis,
      2. verify the split quartic solves,
      3. build the coordinate map
           zeta + {zeta,W3} + 1/2 {{zeta,W3},W3} + {zeta,W4Total},
      4. then inspect single-frequency, deep-water, and toy-broadband outputs.

    This observable map is equivalent to the successful 15-script branch
    {-1, 1, -1}; the old overall sign on W4Back is absorbed into the new
    definition of W4Total, which is required to satisfy

      {H2, W4Total} = -S4RawNonres.
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
    "resonantTerms" -> resonant,
    "nonresTerms" -> nonres,
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

sigmaSimplify[expr_] := FullSimplify[FunctionExpand[expr /. rulesSigma]];

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

broadbandRules = {
  z[1] -> 9 eps/40,   z[-1] -> 9 eps/40,
  z[2] -> 3 eps/40,   z[-2] -> 3 eps/40,
  z[3] -> eps/20,     z[-3] -> eps/20,
  p[1] -> -I (9 eps/40)/Sqrt[T[1]],  p[-1] -> I (9 eps/40)/Sqrt[T[1]],
  p[2] -> -I (3 eps/40)/Sqrt[T[2]],  p[-2] -> I (3 eps/40)/Sqrt[T[2]],
  p[3] -> -I (eps/20)/Sqrt[T[3]],    p[-3] -> I (eps/20)/Sqrt[T[3]]
};

numericStateRules = Flatten@Table[
   {
     aa[k] -> N[(1 + 0.3 k)/(10 + k^2), 30],
     bb[k] -> N[(1 - 0.2 k)/(12 + k^2), 30]
   },
   {k, modes}
];

sigmaSamples = {sigma -> 1/5, sigma -> 2/3, sigma -> 9/10};
numericSigmaRule = sigma -> N[Tanh[1.0], 30];

exprToCubic[expr_, rules_] := Expand[Coefficient[Expand[expr /. rules], eps, 3]];

h3OnlyTarget =
  (3 (1 + sigma^2) (54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 - 4 sigma^8 - 3 sigma^10))/
    (64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));
stokesC33 = (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6);

S4H4Split = splitResonance[Expand[S4H4 /. toNormalRules]];
S4W3Split = splitResonance[Expand[S4W3Induced /. toNormalRules]];
S4RawSplit = splitResonance[Expand[S4Raw /. toNormalRules]];

S4H4Nonres = Expand[S4H4Split["nonres"]];
S4W3Nonres = Expand[S4W3Split["nonres"]];
S4RawNonres = Expand[S4RawSplit["nonres"]];

W4H4 = Expand[solveHomologicalForNegativeSource[S4H4Nonres] /. fromNormalRules];
W4W3 = Expand[solveHomologicalForNegativeSource[S4W3Nonres] /. fromNormalRules];
W4Total = Expand[W4H4 + W4W3];

residualW3Cubic = sigmaSimplify[Expand[H3 - bracketRaw[H2, W3]]];
residualW4H4 = Expand[(bracketRaw[H2, W4H4] + S4H4Nonres /. fromNormalRules) /. toNormalRules];
residualW4W3 = Expand[(bracketRaw[H2, W4W3] + S4W3Nonres /. fromNormalRules) /. toNormalRules];
residualW4Total = Expand[(bracketRaw[H2, W4Total] + S4RawNonres /. fromNormalRules) /. toNormalRules];

residualW4H4Terms = coefficientCheck[residualW4H4];
residualW4W3Terms = coefficientCheck[residualW4W3];
residualW4TotalTerms = coefficientCheck[residualW4Total];

residualW4H4Samples = Table[Quiet[N[residualW4H4 /. rulesSigma /. s /. numericStateRules, 30]], {s, sigmaSamples}];
residualW4W3Samples = Table[Quiet[N[residualW4W3 /. rulesSigma /. s /. numericStateRules, 30]], {s, sigmaSamples}];
residualW4TotalSamples = Table[Quiet[N[residualW4Total /. rulesSigma /. s /. numericStateRules, 30]], {s, sigmaSamples}];

z3W3First = Expand[D[W3, p[-3]]];
z3W3Second = Expand[bracketRaw[z3W3First, W3]];
z3Map = z[3] + z3W3First + 1/2 z3W3Second + Expand[D[W4Total, p[-3]]];

singleC33 =
  sigmaSimplify[2 Coefficient[exprToCubic[z3Map, singleFrequencyRules], eps, 0]];
h3OnlyC33 =
  sigmaSimplify[2 Coefficient[exprToCubic[z[3] + z3W3First + 1/2 z3W3Second, singleFrequencyRules], eps, 0]];

toyObservable = Quiet[N[exprToCubic[z3Map, broadbandRules] /. rulesSigma /. numericSigmaRule, 30]];
toyHamiltonianResidual =
  Quiet[N[exprToCubic[residualW4Total /. fromNormalRules, broadbandRules] /. numericSigmaRule, 30]];

Print["Quartic validation from frozen W3"];
Print[""];
Print["1. W3 cubic Hamiltonian check"];
Print["  H3 - {H2,W3} = ", residualW3Cubic];
Print[""];
Print["2. W4H4 residual"];
Print["  sampled residuals = ", residualW4H4Samples];
Print["  exact residual term count = ", Length[residualW4H4Terms]];
Print["  sample terms = ", Take[residualW4H4Terms, UpTo[6]]];
Print[""];
Print["3. W4W3 residual"];
Print["  sampled residuals = ", residualW4W3Samples];
Print["  exact residual term count = ", Length[residualW4W3Terms]];
Print["  sample terms = ", Take[residualW4W3Terms, UpTo[6]]];
Print[""];
Print["4. W4Total residual"];
Print["  sampled residuals = ", residualW4TotalSamples];
Print["  exact residual term count = ", Length[residualW4TotalTerms]];
Print["  sample terms = ", Take[residualW4TotalTerms, UpTo[6]]];
Print[""];
Print["5. Coordinate map"];
Print["  zetaMap = zeta + {zeta,W3} + 1/2 {{zeta,W3},W3} + {zeta,W4Total}"];
Print["  H3-only C33 from frozen W3 = ", h3OnlyC33];
Print["  H3-only target             = ", h3OnlyTarget];
Print["  H3-only residual           = ", sigmaSimplify[h3OnlyC33 - h3OnlyTarget]];
Print[""];
Print["6. Single-frequency C33"];
Print["  full quartic C33 = ", singleC33];
Print["  Stokes/VWA C33   = ", stokesC33];
Print["  residual         = ", sigmaSimplify[singleC33 - stokesC33]];
Print[""];
Print["7. Deep-water limit"];
Print["  sigma -> 1 gives ", FullSimplify[singleC33 /. sigma -> 1]];
Print[""];
Print["8. Toy broadband state"];
Print["  z[3] cubic observable on toy state = ", toyObservable];
Print["  Hamiltonian residual on same state = ", toyHamiltonianResidual];
Print[""];
Print["Interpretation:"];
Print["  If the Hamiltonian residual above is zero but the toy broadband"];
Print["  observable still looks unsatisfactory, the remaining issue is"];
Print["  observable content, truncation, or higher-order structure rather"];
Print["  than the quartic homological solve itself."];
