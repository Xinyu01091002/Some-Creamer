(* ::Package:: *)
(*  Iterative finite-depth H3/H4 transform search.

    Strategy
    --------
    We keep the accepted direct-H4 / non-resonant-normal-form structure,
    but search the remaining sign conventions step by step:

      1. pick the W3 observable sign that reproduces the known H3-only C33;
      2. build K4 = H4 + bracketSign 1/2 {H3,W3};
      3. absorb only the non-resonant quartic part through W4;
      4. check the quartic homological residual coefficient-by-coefficient;
      5. among the surviving candidates, check the single-frequency C33;
      6. print toy-broadband diagnostics for the surviving transform(s).

    This avoids one giant FullSimplify on the whole quartic polynomial and
    follows the Creamer-style logic more closely: first verify the
    Hamiltonian identity, then verify the induced coordinate map.
*)

ClearAll["Global`*"];

$Assumptions = (And @@ Table[T[n] > 0, {n, 1, 6}]) &&
  Element[sigma, Reals] && 0 < sigma < 1;

modes = {-3, -2, -1, 1, 2, 3};

theta[0] := 0;
theta[n_Integer] /; n != 0 := T[Abs[n]];
omega[n_Integer] /; n != 0 := Sqrt[theta[n]];
ksq[n_Integer] := n^2;

logLine[msg_] := (WriteString[$Output, msg <> "\n"]; Flush[$Output];);

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
    "resonantTerms" -> resonant,
    "nonresTerms" -> nonres,
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

broadbandRules = {
  z[1] -> 9 eps/40,   z[-1] -> 9 eps/40,
  z[2] -> 3 eps/40,   z[-2] -> 3 eps/40,
  z[3] -> eps/20,     z[-3] -> eps/20,
  p[1] -> -I (9 eps/40)/Sqrt[T[1]],  p[-1] -> I (9 eps/40)/Sqrt[T[1]],
  p[2] -> -I (3 eps/40)/Sqrt[T[2]],  p[-2] -> I (3 eps/40)/Sqrt[T[2]],
  p[3] -> -I (eps/20)/Sqrt[T[3]],    p[-3] -> I (eps/20)/Sqrt[T[3]]
};

numericSigmaRule = sigma -> N[Tanh[1.0], 30];
numericSigmaSamples = {sigma -> 1/5, sigma -> 2/3, sigma -> 9/10};
stokesC33 = (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6);
h3OnlyTarget =
  (3 (1 + sigma^2) (54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 - 4 sigma^8 - 3 sigma^10))/
    (64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

exprToCubic[expr_, rules_] := Expand[Coefficient[Expand[expr /. rules], eps, 3]];

monomialFromPowers[powers_] := Times @@ MapThread[#1^#2 &, {normalVars, powers}];

coefficientCheck[expr_] := Module[{rules, numericFiltered, checked},
  rules = CoefficientRules[Expand[expr], normalVars];
  numericFiltered = Select[
    rules,
    Module[{vals},
      vals = N[Table[#[[2]] /. rulesSigma /. s, {s, numericSigmaSamples}], 40];
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

logLine["14: checking cubic Hamiltonian sign."];
r3Data = Table[
   <|
     "hs3" -> hs3,
     "residualTerms" -> coefficientCheck[H3 + hs3 bracketRaw[H2, W3]]
   |>,
   {hs3, {-1, 1}}
];
r3Good = Select[r3Data, Length[#["residualTerms"]] == 0 &];

z3W3First = Expand[D[W3, p[-3]]];
z3W3Second = Expand[bracketRaw[z3W3First, W3]];
z3H3Map[s3_] := z[3] + s3 z3W3First + 1/2 s3^2 z3W3Second;

logLine["14: checking W3 observable sign against H3-only C33."];
s3Data = Table[
   <|
     "s3" -> s3,
     "c33" -> Simplify[2 Coefficient[exprToCubic[z3H3Map[s3], singleFrequencyRules] /. rulesSigma, eps, 0]],
     "residual" -> Simplify[
       2 Coefficient[exprToCubic[z3H3Map[s3], singleFrequencyRules] /. rulesSigma, eps, 0] - h3OnlyTarget]
   |>,
   {s3, {-1, 1}}
];
s3Good = Select[s3Data, TrueQ[#["residual"] === 0] &];
chosenS3 = If[Length[s3Good] > 0, First[s3Good]["s3"], 1];

searchSeeds = {
  <|"bracketSign" -> -1, "inverseSign" -> -1, "s4" -> 1|>,
  <|"bracketSign" -> -1, "inverseSign" -> 1, "s4" -> -1|>
};

logLine["14: entering quartic candidate scan."];
candidateData = Table[
   Module[
     {seed, k4Orig, k4Normal, split, w4Normal, w4Back, adH2W4, r4Expr,
      r4Terms, z3Map, c33, c33Residual, toyObservable, toyHamResidual},
     seed = searchSeeds[[idx]];
     logLine[
      "14: candidate " <> ToString[idx] <> "/" <> ToString[Length[searchSeeds]] <>
      " " <> ToString[{seed["bracketSign"], seed["inverseSign"], seed["s4"]}, InputForm]
     ];
     logLine["14:   building K4."];
     k4Orig = Expand[H4Direct + seed["bracketSign"] 1/2 bracketRaw[H3, W3]];
     k4Normal = Expand[k4Orig /. toNormalRules];
     logLine["14:   splitting resonant/non-resonant terms."];
     split = splitResonance[k4Normal];
     logLine["14:   solving W4 from non-resonant quartic part."];
     w4Normal = seed["inverseSign"] solveHomologicalNonres[split["nonres"]];
     logLine["14:   mapping W4 back to (z,p)."];
     w4Back = Expand[w4Normal /. fromNormalRules];
     logLine["14:   checking quartic homological residual."];
     adH2W4 = Expand[(bracketRaw[H2, w4Back]) /. toNormalRules];
     r4Expr = Expand[split["nonres"] - adH2W4];
     r4Terms = coefficientCheck[r4Expr];
     logLine["14:   checking observable-level single-frequency and toy broadband diagnostics."];
     z3Map = z3H3Map[chosenS3] + seed["s4"] Expand[D[w4Back, p[-3]]];
     c33 = Simplify[2 Coefficient[exprToCubic[z3Map, singleFrequencyRules] /. rulesSigma, eps, 0]];
     c33Residual = Simplify[c33 - stokesC33];
     toyObservable = N[exprToCubic[z3Map, broadbandRules] /. rulesSigma /. numericSigmaRule, 30];
     toyHamResidual = N[exprToCubic[(r4Expr /. fromNormalRules), broadbandRules] /. numericSigmaRule, 30];
     logLine["14:   candidate complete."];
     <|
       "bracketSign" -> seed["bracketSign"],
       "inverseSign" -> seed["inverseSign"],
       "s4" -> seed["s4"],
       "K4ResCount" -> split["resonantCount"],
       "K4NonresCount" -> split["nonresCount"],
       "R4ResidualCount" -> Length[r4Terms],
       "R4ResidualTerms" -> Take[r4Terms, UpTo[8]],
       "singleC33" -> c33,
       "singleResidual" -> c33Residual,
       "toyObservable" -> toyObservable,
       "toyHamiltonianResidual" -> toyHamResidual
     |>
   ],
   {idx, Length[searchSeeds]}
];

goodCandidates = Select[
  candidateData,
  Length[#["R4ResidualTerms"]] == 0 && TrueQ[#["singleResidual"] === 0] &
];

Print["Iterative finite-depth H3/H4 transform search"];
Print[""];
Print["Cubic Hamiltonian sign candidates hs3 from H3 + hs3 {H2,W3}: "];
Print[r3Data];
Print["Exact cubic Hamiltonian candidates: ", r3Good];
Print[""];

Print["Observable sign scan for W3 using H3-only single-frequency C33:"];
Print[s3Data];
Print["Chosen s3: ", chosenS3];
Print[""];

Print["Quartic candidate scan (absorb only non-resonant K4 terms):"];
Do[
  Print[
    "  {bracketSign,inverseSign,s4} = ",
    {row["bracketSign"], row["inverseSign"], row["s4"]},
    " -> K4{res,nonres} = ",
    {row["K4ResCount"], row["K4NonresCount"]},
    "; R4 residual count = ", row["R4ResidualCount"],
    "; single residual = ", row["singleResidual"],
    "; toy observable = ", row["toyObservable"],
    "; toy Hamiltonian residual = ", row["toyHamiltonianResidual"]
  ],
  {row, candidateData}
];
Print[""];

Print["Candidates satisfying both quartic homological cancellation and single-frequency C33:"];
Print[goodCandidates];
Print[""];

If[Length[goodCandidates] > 0,
  Print["Sample residual terms for surviving candidates (should be empty):"];
  Do[
    Print[
      "  candidate ",
      {row["bracketSign"], row["inverseSign"], row["s4"]},
      " residual terms: ",
      row["R4ResidualTerms"]
    ],
    {row, goodCandidates}
  ],
  Print["No exact surviving candidate found in this sign family."];
  Print["Sample residual terms from the best candidate by residual count:"];
  With[{best = First@SortBy[candidateData, {#["R4ResidualCount"] &, ToString[#["singleResidual"], InputForm] &}]},
    Print["  best candidate = ", {best["bracketSign"], best["inverseSign"], best["s4"]}];
    Print["  first residual terms = ", best["R4ResidualTerms"]];
  ]
];
