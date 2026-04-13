(* ::Package:: *)
(*  Freeze the accepted finite-depth W3 basis before entering quartic work.

    This script locks two sign facts already diagnosed in 10/13/15:

      1. Hamiltonian convention:
           H3 - {H2, W3} = 0

      2. Observable convention:
           zeta + {zeta, W3} + 1/2 {{zeta, W3}, W3}
         reproduces the current H3-only single-frequency C33 target.

    Once these checks pass, later quartic scripts treat this W3 as frozen.
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

exprToCubic[expr_, rules_] := Expand[Coefficient[Expand[expr /. rules], eps, 3]];

h3OnlyTarget =
  (3 (1 + sigma^2) (54 + 39 sigma^2 + 70 sigma^4 - 28 sigma^6 - 4 sigma^8 - 3 sigma^10))/
    (64 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

r3Data = Table[
   <|
     "hs3" -> hs3,
     "residual" -> sigmaSimplify[Expand[H3 + hs3 bracketRaw[H2, W3]]]
   |>,
   {hs3, {-1, 1}}
];

z3W3First = Expand[D[W3, p[-3]]];
z3W3Second = Expand[bracketRaw[z3W3First, W3]];
z3H3Map = z[3] + z3W3First + 1/2 z3W3Second;
z3H3MapC33 = sigmaSimplify[2 Coefficient[exprToCubic[z3H3Map, singleFrequencyRules], eps, 0]];

Print["Frozen finite-depth W3 basis check"];
Print[""];
Print["Hamiltonian sign scan for H3 + hs3 {H2,W3}:"]; 
Do[
  Print["  hs3 = ", row["hs3"], " -> residual = ", row["residual"]],
  {row, r3Data}
];
Print["  frozen Hamiltonian sign: hs3 = -1"];
Print["  exact cubic check H3 - {H2,W3} == 0: ",
  TrueQ[sigmaSimplify[Expand[H3 - bracketRaw[H2, W3]]] === 0]
];
Print[""];
Print["Observable sign lock for the frozen W3 basis:"];
Print["  s3 = +1"];
Print["  z[3] map = z[3] + {z[3],W3} + 1/2 {{z[3],W3},W3}"];
Print["  H3-only C33 from frozen W3 = ", z3H3MapC33];
Print["  current H3-only target      = ", h3OnlyTarget];
Print["  residual                    = ", sigmaSimplify[z3H3MapC33 - h3OnlyTarget]];
