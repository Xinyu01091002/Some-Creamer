(* ::Package:: *)
(*
  Debug the relation between the 01 lambda-flow quadratic vector field and a
  W3 bracket representation.

  The previous Hamiltonian-normal-form coordinate-map check used a bracket
  generator W3, but it did not reproduce the H3-only C33 coefficient from
  01_single_frequency_limit.wl.  This script compares the two constructions
  directly at the vector-field level before adding any H4/cubic terms.
*)

ClearAll["Global`*"];

$Assumptions = (And @@ Table[T[n] > 0, {n, 1, 6}]) &&
  Element[sigma, Reals] && 0 < sigma < 1;

modes = {-3, -2, -1, 1, 2, 3};

theta[0] := 0;
theta[n_Integer] /; n != 0 := T[Abs[n]];
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

Dker[dest_Integer, p_Integer, r_Integer] :=
  Dker[dest, p, r] =
    dFinite[theta[dest], theta[p], theta[r], ksq[dest], ksq[p], ksq[r]];

BzKer[dest_Integer, p_Integer, r_Integer] :=
  BzKer[dest, p, r] =
    bFinite[theta[p], theta[r], theta[dest], ksq[p], ksq[r], ksq[dest]];

BphiKer[dest_Integer, p_Integer, r_Integer] :=
  BphiKer[dest, p, r] =
    bFinite[theta[dest], theta[p], theta[r], ksq[dest], ksq[p], ksq[r]];

bracketRaw[F_, G_] :=
  Sum[
    D[F, z[k]] D[G, p[-k]] - D[F, p[k]] D[G, z[-k]],
    {k, modes}
  ];

(* Direct 01 lambda-flow vector field, as a derivation on arbitrary F. *)
lambdaDz[n_] := Sum[
  If[p1 + r == n,
    3 Dker[n, p1, r] p[p1] p[r] + BzKer[n, p1, r] z[p1] z[r],
    0
  ],
  {p1, modes}, {r, modes}
];

lambdaDp[n_] := Sum[
  If[p1 + r == n,
    -2 BphiKer[n, p1, r] z[p1] p[r],
    0
  ],
  {p1, modes}, {r, modes}
];

lambdaDeriv[F_] := Sum[
  D[F, z[n]] lambdaDz[n] + D[F, p[n]] lambdaDp[n],
  {n, modes}
];

(* Candidate generator used in 09. *)
W3Candidate = Sum[
  If[a + b + c == 0,
    Dker[a, b, c] p[a] p[b] p[c]
    + BphiKer[a, b, c] z[a] z[b] p[c],
    0
  ],
  {a, modes}, {b, modes}, {c, modes}
];

candidateDz[n_] := bracketRaw[z[n], W3Candidate];
candidateDp[n_] := bracketRaw[p[n], W3Candidate];

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
  T[3] -> 3 Tanh[3 ArcTanh[sigma]]
};

cosC33FromPositiveMode[expr_] :=
  FullSimplify[2 Coefficient[Expand[expr /. linearWaveRules], eps, 3] /. rulesSigma];

lambdaSecondOrderZ2 = lambdaDeriv[z[2]];
lambdaThirdOrderZ3ViaFlow = 1/2 lambdaDeriv[lambdaDeriv[z[3]]];

candidateSecondOrderZ2 = candidateDz[2];
candidateThirdOrderZ3ViaFlow = 1/2 bracketRaw[bracketRaw[z[3], W3Candidate], W3Candidate];

creamerH3OnlyC33 =
  (3 (1 + sigma^2) (54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 -
      4 sigma^8 - 3 sigma^10))/
    (64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

lambdaC33 = cosC33FromPositiveMode[lambdaThirdOrderZ3ViaFlow];
candidateC33 = cosC33FromPositiveMode[candidateThirdOrderZ3ViaFlow];

Print["Lambda-flow W3 alignment debug"];
Print[""];
Print["Compare vector fields on z[2] after linear single-frequency substitution:"];
Print["  direct lambda Dz[2]:     ",
  FullSimplify[lambdaSecondOrderZ2 /. linearWaveRules /. rulesSigma]
];
Print["  bracket candidate Dz[2]: ",
  FullSimplify[candidateSecondOrderZ2 /. linearWaveRules /. rulesSigma]
];
Print[""];
Print["C33 from applying the quadratic vector field twice to z[3]:"];
Print["  direct lambda-flow:      ", lambdaC33];
Print["  bracket candidate:       ", candidateC33];
Print["  01 H3-only target:       ", creamerH3OnlyC33];
Print["  lambda - target:         ", FullSimplify[lambdaC33 - creamerH3OnlyC33]];
Print["  candidate - target:      ", FullSimplify[candidateC33 - creamerH3OnlyC33]];
Print[""];
Print["Raw mismatch in vector field for modes +/-2, +/-3:"];
Do[
  Print[
    "  n=", n,
    " Dz diff=", FullSimplify[(candidateDz[n] - lambdaDz[n]) /. linearWaveRules],
    " Dp diff=", FullSimplify[(candidateDp[n] - lambdaDp[n]) /. linearWaveRules]
  ],
  {n, {-3, -2, 2, 3}}
];
