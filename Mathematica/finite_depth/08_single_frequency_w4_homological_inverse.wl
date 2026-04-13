(* ::Package:: *)
(*
  Single-frequency W4 homological-inverse check.

  This is the next step after 07_single_frequency_quartic_normal_form_attempt:

      K4 = H4 + 1/2 {H3, W3}
      {H2, W4} = K4_nonres

  Instead of solving the full four-wave kernel, this script diagonalizes H2
  in the truncated 1-D mode set +/-1,+/-2,+/-3, solves the homological equation
  monomial-by-monomial in linear oscillator variables, maps W4 back to the
  original (zeta, phi_s) Fourier variables, and extracts the coefficient of
  z[1]^3 p[-3].  That is the minimal coefficient which would change the
  single-frequency third harmonic z[3] through {z[3], W4}.

  Conventions:
    bracketRaw[F,G] is {F,G}.
    a_k = z_k + I omega_k p_k, b_k = z_k - I omega_k p_k.
    {H2, M} = I (sum_a omega - sum_b omega) M for a normal monomial M.

  The sign of the final surface increment is therefore still a convention
  check against the lambda-flow scripts; the script prints the main sign
  variants rather than silently choosing one.
*)

ClearAll["Global`*"];

$Assumptions = Element[sigma, Reals] && 0 < sigma < 1;

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

Bker[a_, b_, c_] := Bker[a, b, c] =
  bFinite[theta[a], theta[b], theta[c], ksq[a], ksq[b], ksq[c]];

H3ker[a_, b_, c_] :=
  1/4 (ksq[b] + ksq[c] - ksq[a] - 2 theta[b] theta[c]);

H4OrderedKer[a_, b_, c_, d_] :=
  1/2 theta[a] theta[c + d] theta[d]
  - 1/4 (theta[a] ksq[d] + ksq[a] theta[d]);

bracketRaw[F_, G_] :=
  Sum[
    D[F, z[k]] D[G, p[-k]] - D[F, p[k]] D[G, z[-k]],
    {k, modes}
  ];

H2 = 1/2 Sum[
    z[k] z[-k] + theta[k] p[k] p[-k],
    {k, modes}
  ];

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
      + Bker[a, b, c] z[a] z[b] p[c],
      0
    ],
    {a, modes}, {b, modes}, {c, modes}
  ];

H4 = Sum[
    If[a + b + c + d == 0,
      z[b] z[c] p[a] p[d] H4OrderedKer[a, b, c, d],
      0
    ],
    {a, modes}, {b, modes}, {c, modes}, {d, modes}
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

solveHomological[poly_] := Module[{terms, mon, coeff, weight, solvedTerms},
  terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
  solvedTerms = Table[
    mon = Times @@ Table[var^Exponent[term, var], {var, normalVars}];
    coeff = Together[term/mon];
    weight = normalWeight[mon];
    If[TrueQ[FullSimplify[weight == 0]],
      0,
      coeff mon/(I weight)
    ],
    {term, terms}
  ];
  Total[solvedTerms]
];

computeW4Coeff[bracketSign_, w3Sign_, inverseSign_] := Module[
  {K4Local, K4NormalLocal, W4NormalLocal, W4BackLocal},
  K4Local = Expand[H4 + bracketSign 1/2 bracketRaw[H3, w3Sign W3]];
  K4NormalLocal = Expand[K4Local /. toNormalRules];
  W4NormalLocal = inverseSign solveHomological[K4NormalLocal];
  W4BackLocal = Expand[W4NormalLocal /. fromNormalRules];
  Together[Coefficient[W4BackLocal, z[1]^3 p[-3]]]
];

variantDataRaw = Flatten[
  Table[
    <|
      "bracketSign" -> bs,
      "w3Sign" -> ws,
      "inverseSign" -> is,
      "coeff" -> computeW4Coeff[bs, ws, is]
    |>,
    {bs, {-1, 1}}, {ws, {-1, 1}}, {is, {-1, 1}}
  ],
  2
];

w4CoeffRaw = computeW4Coeff[1, 1, 1];

rulesSigma = {
  T[1] -> sigma,
  T[2] -> 2 Tanh[2 ArcTanh[sigma]],
  T[3] -> 3 Tanh[3 ArcTanh[sigma]],
  T[4] -> 4 Tanh[4 ArcTanh[sigma]],
  T[5] -> 5 Tanh[5 ArcTanh[sigma]],
  T[6] -> 6 Tanh[6 ArcTanh[sigma]]
};

w4CoeffSigma = FullSimplify[w4CoeffRaw /. rulesSigma];
w4CoeffDeepWater = FullSimplify[w4CoeffSigma /. sigma -> 1];
w4FiniteDepthOnlySigma = FullSimplify[w4CoeffSigma - w4CoeffDeepWater];

triadCheckMonomial = z[-2] p[1] p[1];
triadH3Coeff = Together[Coefficient[H3, triadCheckMonomial]];
triadAdCoeff = Together[Coefficient[bracketRaw[H2, W3], triadCheckMonomial]];
triadRatio = FullSimplify[triadAdCoeff/triadH3Coeff /. rulesSigma];

variantDataSigma = Map[
  Append[#,
    "coeffSigma" -> FullSimplify[#["coeff"] /. rulesSigma]
  ] &,
  variantDataRaw
];

variantDataSummary = Map[
  Append[#,
    "deepWater" -> FullSimplify[#["coeffSigma"] /. sigma -> 1]
  ] &,
  variantDataSigma
];

requiredW4SurfaceCoeff =
  (3 (-1 + sigma^2)^2 (27 + 60 sigma^2 + 50 sigma^4 + 4 sigma^6 + 3 sigma^8))/
    (16 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

residualPlus = FullSimplify[w4CoeffSigma - requiredW4SurfaceCoeff];
residualFiniteDepthOnly = FullSimplify[w4FiniteDepthOnlySigma - requiredW4SurfaceCoeff];

variantResidualSummary = Map[
  Append[#,
    "residualSigma" -> FullSimplify[#["coeffSigma"] - requiredW4SurfaceCoeff]
  ] &,
  variantDataSummary
];

Print["Single-frequency W4 homological-inverse check"];
Print[""];
Print["Mode truncation: +/-1,+/-2,+/-3"];
Print["Normal variables: a_k = z_k + I omega_k p_k, b_k = z_k - I omega_k p_k"];
Print["Default variant solves {H2,W4}=K4_nonres with K4=H4+1/2 {H3,W3}."];
Print[""];
Print["Cubic generator calibration check on monomial z[-2] p[1]^2:"];
Print["  coeff in H3:       ", triadH3Coeff];
Print["  coeff in {H2,W3}:  ", triadAdCoeff];
Print["  ratio in sigma:    ", triadRatio];
Print[""];
Print["Coefficient of z[1]^3 p[-3] in W4 after mapping back:"];
Print["  raw T[n] form: ", w4CoeffRaw];
Print["  sigma form:    ", w4CoeffSigma];
Print["  finite-depth-only part W4(sigma)-W4(1): ", w4FiniteDepthOnlySigma];
Print[""];
Print["Required calibrated coefficient from finite-depth Stokes C[3,3]:"];
Print["  ", requiredW4SurfaceCoeff];
Print[""];
Print["Residuals against required coefficient:"];
Print["  W4 - required:  ", residualPlus];
Print["  (W4-W4_deep_water) - required: ", residualFiniteDepthOnly];
Print[""];
Print["Deep-water checks sigma -> 1:"];
Print["  {default W4, required}: ",
  FullSimplify[{w4CoeffSigma, requiredW4SurfaceCoeff} /. sigma -> 1]
];
Print[""];
Print["Sign-convention scan:"];
Print["  bracketSign is the sign in H4 + bracketSign/2 {H3, w3Sign W3}."];
Print["  inverseSign is the sign multiplying the homological inverse."];
Do[
  Print[
    "  {bracketSign,w3Sign,inverseSign} = ",
    {row["bracketSign"], row["w3Sign"], row["inverseSign"]},
    " -> deep water coeff ",
    row["deepWater"],
    "; coeff - required = ",
    row["residualSigma"]
  ],
  {row, variantResidualSummary}
];
