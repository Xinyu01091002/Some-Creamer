(* ::Package:: *)
(*  Overall finite-depth H3/H4 homological-consistency check.

    Purpose
    -------
    The single-frequency C33 checks in 11/12 show that direct H4 plus a
    non-resonant W4 can recover the finite-depth Stokes/VWA third-harmonic
    coefficient.  This script asks the more structural question:

      - does the chosen W3/W4 convention satisfy the cubic/quartic
        homological equations as Hamiltonian identities?
      - once that convention is fixed, what does the same map predict for a
        small broadband state on the truncated mode set +/-1,+/-2,+/-3?

    The goal is to separate:
      (a) a wrong sign/gauge/convention,
      (b) a failure of the non-resonant homological solve itself,
      (c) a genuinely broadband observable mismatch not diagnosed by the
          single-frequency C33 test.
*)

ClearAll["Global`*"];

$Assumptions = (And @@ Table[T[n] > 0, {n, 1, 6}]) &&
  Element[sigma, Reals] && 0 < sigma < 1;

modes = {-3, -2, -1, 1, 2, 3};

theta[0] := 0;
theta[n_Integer] /; n != 0 := T[Abs[n]];
omega[n_Integer] /; n != 0 := Sqrt[theta[n]];
gamma[n_Integer] /; n != 0 := Sqrt[theta[n]];
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
    If[TrueQ[FullSimplify[weight == 0]],
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
    If[TrueQ[FullSimplify[weight == 0]], 0, coeff mon/(I weight)],
    {term, terms}
  ];
  Total[solved]
];

polySummary[expr_, n_: 8] := Module[{terms},
  terms = If[TrueQ[Expand[expr] === 0], {}, List @@ Expand[expr]];
  <|
    "count" -> Length[terms],
    "sample" -> Take[terms, UpTo[n]]
  |>
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

(* Small symmetric broadband toy state on +/-1,+/-2,+/-3. *)
broadbandRules = {
  z[1] -> 9 eps/40,   z[-1] -> 9 eps/40,
  z[2] -> 3 eps/40,   z[-2] -> 3 eps/40,
  z[3] -> eps/20,     z[-3] -> eps/20,
  p[1] -> -I (9 eps/40)/Sqrt[T[1]],  p[-1] -> I (9 eps/40)/Sqrt[T[1]],
  p[2] -> -I (3 eps/40)/Sqrt[T[2]],  p[-2] -> I (3 eps/40)/Sqrt[T[2]],
  p[3] -> -I (eps/20)/Sqrt[T[3]],    p[-3] -> I (eps/20)/Sqrt[T[3]]
};

exprToCubic[expr_, rules_] := Expand[Coefficient[Expand[expr /. rules], eps, 3]];

numericSigmaRule = sigma -> N[Tanh[1.0], 30];

logLine[msg_] := (WriteString[$Output, msg <> "\n"]; Flush[$Output];);

(* --- cubic sign scan --- *)
r3Data = Table[
   With[{res = Simplify[Expand[H3 + hs3 bracketRaw[H2, W3]] /. rulesSigma]},
     <|
       "hs3" -> hs3,
       "residual" -> res,
       "summary" -> polySummary[res]
     |>
   ],
   {hs3, {-1, 1}}
];
r3Zero = Select[r3Data, TrueQ[#["residual"] === 0] &];

z3W3First = Expand[D[W3, p[-3]]];
z3W3Second = Expand[bracketRaw[z3W3First, W3]];

z3H3Map[s3_] :=
  z[3] + s3 z3W3First + 1/2 s3^2 z3W3Second;

h3OnlyTarget =
  (3 (1 + sigma^2) (54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 - 4 sigma^8 - 3 sigma^10))/
    (64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

s3MapData = Table[
   <|
     "s3" -> s3,
     "c33" -> FullSimplify[2 Coefficient[exprToCubic[z3H3Map[s3], singleFrequencyRules] /. rulesSigma, eps, 0]],
     "residualToH3Only" -> Simplify[
       2 Coefficient[exprToCubic[z3H3Map[s3], singleFrequencyRules] /. rulesSigma, eps, 0] - h3OnlyTarget]
   |>,
   {s3, {-1, 1}}
];

chosenS3 =
  If[Length[Select[s3MapData, TrueQ[#["residualToH3Only"] === 0] &]] > 0,
    First[Select[s3MapData, TrueQ[#["residualToH3Only"] === 0] &]]["s3"],
    1
  ];

(* --- quartic scan --- *)
quarticConventionSeeds = {
  <|"bracketSign" -> -1, "inverseSign" -> -1, "s4" -> 1|>,
  <|"bracketSign" -> -1, "inverseSign" -> 1, "s4" -> -1|>
};

logLine["13: finished cubic setup; entering quartic base solves."];
quarticBaseData = Flatten[
   Table[
     Module[
       {k4Orig, k4Normal, split, w4Normal, w4Back, adH2W4Normal, r4Natural},
       logLine[
        "13: solving quartic base for {bracketSign,inverseSign,s4} = " <>
        ToString[{seed["bracketSign"], seed["inverseSign"], seed["s4"]}, InputForm]
       ];
       k4Orig = Expand[H4Direct + seed["bracketSign"] 1/2 bracketRaw[H3, W3]];
       k4Normal = Expand[k4Orig /. toNormalRules];
       split = splitResonance[k4Normal];
       w4Normal = seed["inverseSign"] solveHomologicalNonres[k4Normal];
       w4Back = Expand[w4Normal /. fromNormalRules];
       adH2W4Normal = Expand[(bracketRaw[H2, w4Back]) /. toNormalRules];
       r4Natural = Simplify[Expand[split["nonres"] - adH2W4Normal] /. rulesSigma];
       logLine[
        "13: finished quartic base for {bracketSign,inverseSign,s4} = " <>
        ToString[{seed["bracketSign"], seed["inverseSign"], seed["s4"]}, InputForm]
       ];
       <|
         "bracketSign" -> seed["bracketSign"],
         "inverseSign" -> seed["inverseSign"],
         "s4Preferred" -> seed["s4"],
         "K4AllCount" -> polySummary[k4Normal]["count"],
         "K4ResCount" -> split["resonantCount"],
         "K4NonresCount" -> split["nonresCount"],
         "K4Normal" -> k4Normal,
         "split" -> split,
         "W4Back" -> w4Back,
         "R4Natural" -> r4Natural,
         "R4Summary" -> polySummary[r4Natural]
       |>
     ],
     {seed, quarticConventionSeeds}
   ],
   1
];

quarticData = Flatten[
   Table[
     Module[{z3W4First, z3Overall, singleResidual, toyValue, toyResidualNorm},
       z3W4First = Expand[D[row["W4Back"], p[-3]]];
       z3Overall = z3H3Map[chosenS3] + row["s4Preferred"] z3W4First;
       singleResidual = Simplify[
         2 Coefficient[exprToCubic[z3Overall, singleFrequencyRules] /. rulesSigma, eps, 0]
         - (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6)
       ];
       toyValue = N[(exprToCubic[z3Overall, broadbandRules] /. rulesSigma /. numericSigmaRule), 30];
       toyResidualNorm = N[(exprToCubic[row["R4Natural"] /. fromNormalRules, broadbandRules] /. numericSigmaRule), 30];
       Association[
         KeyDrop[row, {"K4Normal", "split", "W4Back"}],
         <|
           "s3" -> chosenS3,
           "s4" -> row["s4Preferred"],
           "singleResidual" -> singleResidual,
           "toyZ3Cubic" -> toyValue,
           "toyHamiltonianResidual" -> toyResidualNorm
         |>
       ]
     ],
     {row, quarticBaseData}
   ],
   1
];

bestConventions = Select[
  quarticData,
  TrueQ[#["R4Natural"] === 0] && TrueQ[#["singleResidual"] === 0] &
];

chosen =
  If[Length[bestConventions] > 0,
    First[bestConventions],
    First[SortBy[quarticData, {polySummary[#["R4Natural"]]["count"] &, polySummary[#["singleResidual"]]["count"] &}]]
  ];

chosenBase = First@Select[
  quarticBaseData,
  #["bracketSign"] == chosen["bracketSign"] && #["inverseSign"] == chosen["inverseSign"] &
];
chosenK4Normal = chosenBase["K4Normal"];
chosenSplit = chosenBase["split"];
chosenW4Back = chosenBase["W4Back"];
chosenAdH2W4 = Expand[(bracketRaw[H2, chosenW4Back]) /. toNormalRules];
chosenR4 = Simplify[Expand[chosenSplit["nonres"] - chosenAdH2W4] /. rulesSigma];

z3Chosen = z3H3Map[chosen["s3"]] + chosen["s4"] Expand[D[chosenW4Back, p[-3]]];
singleChosenC33 = Simplify[
  2 Coefficient[exprToCubic[z3Chosen, singleFrequencyRules] /. rulesSigma, eps, 0]
];
broadbandChosenCubic = Expand[exprToCubic[z3Chosen, broadbandRules] /. rulesSigma /. numericSigmaRule];

Print["Finite-depth overall H3/H4 homological-consistency check"];
Print[""];
Print["Cubic Hamiltonian sign scan R3 = H3 + hs3 {H2,W3}:"]; 
Do[
  Print["  hs3 = ", row["hs3"], " -> term count ", row["summary", "count"], "; residual = ", row["residual"]],
  {row, r3Data}
];
Print["  exact zero R3 candidates: ", r3Zero];
Print[""];

Print["Observable-level W3 scan on single-frequency C33:"];
Do[
  Print[
    "  s3 = ", row["s3"],
    " -> C33 = ", row["c33"],
    "; residual to H3-only target = ", row["residualToH3Only"]
  ],
  {row, s3MapData}
];
Print[""];

Print["Quartic convention scan (Hamiltonian residual and single-frequency observable):"];
Do[
  Print[
    "  {bracketSign,inverseSign,s3,s4} = ",
    {row["bracketSign"], row["inverseSign"], row["s3"], row["s4"]},
    " -> K4{all,res,nonres} = ",
    {row["K4AllCount"], row["K4ResCount"], row["K4NonresCount"]},
    "; R4 terms = ", row["R4Summary", "count"],
    "; single-frequency residual = ", row["singleResidual"],
    "; broadband toy z3 cubic = ", row["toyZ3Cubic"],
    "; broadband toy Hamiltonian residual = ", row["toyHamiltonianResidual"]
  ],
  {row, quarticData}
];
Print[""];
Print["Conventions satisfying both R4_natural==0 and single-frequency residual==0:"];
Print[bestConventions];
Print[""];

Print["Chosen convention for downstream diagnostics:"];
Print["  {bracketSign,inverseSign,s3,s4} = ",
  {chosen["bracketSign"], chosen["inverseSign"], chosen["s3"], chosen["s4"]}];
Print["  K4 term counts {all,res,nonres} = ",
  {chosen["K4AllCount"], chosen["K4ResCount"], chosen["K4NonresCount"]}];
Print[""];

Print["Chosen quartic split summaries:"];
Print["  sample resonant terms: ", Take[chosenSplit["resonantTerms"], UpTo[6]]];
Print["  sample non-resonant terms: ", Take[chosenSplit["nonresTerms"], UpTo[6]]];
Print[""];

Print["Chosen R4 = K4_nonres - {H2,W4} (normal variables, finite-depth substitution):"];
Print["  term count: ", polySummary[chosenR4]["count"]];
Print["  residual: ", chosenR4];
Print[""];

Print["Single-frequency observable check under chosen convention:"];
Print["  C33 = ", singleChosenC33];
Print["  residual to Stokes/VWA = ",
  Simplify[singleChosenC33 - (27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6)]];
Print["  deep-water limit = ", Simplify[singleChosenC33 /. sigma -> 1]];
Print[""];

Print["Broadband toy-state observable diagnostics (sigma = tanh(1), small symmetric +/-1,+/-2,+/-3 state):"];
Print["  z3 cubic observable = ", broadbandChosenCubic];
Print["  z3 cubic absolute value = ", N[Abs[broadbandChosenCubic], 20]];
Print["  Hamiltonian residual evaluated on same toy state = ",
  N[(exprToCubic[chosenR4 /. fromNormalRules, broadbandRules] /. numericSigmaRule), 20]];
Print[""];

Print["Interpretation:"];
Print["  If the chosen convention gives R4=0 but broadband toy observables still"];
Print["  look suspicious, the remaining issue is not the quartic homological solve"];
Print["  itself, but the broadband observable content captured by H3/H4 on this"];
Print["  truncated mode set."];
