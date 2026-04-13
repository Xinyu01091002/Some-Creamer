(* ::Package:: *)
(*  Make the quartic source explicit before solving for W4.

    The frozen convention from 16 is:

      H3 - {H2, W3} = 0

    so the quartic source is organized as

      S4Raw = S4H4 + S4W3Induced
            = H4Direct + (-{H3,W3} + 1/2 {{H2,W3},W3}).

    The accepted shorthand H4Direct - 1/2 {H3,W3} is not the starting point
    here; it is recovered only after using the frozen cubic relation.
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

currentShorthand = Expand[S4H4 - 1/2 bracketRaw[H3, W3]];
equivalenceResidual = Expand[S4Raw - currentShorthand];
equivalenceViaFrozenW3 =
  Expand[equivalenceResidual + 1/2 bracketRaw[H3 - bracketRaw[H2, W3], W3]];

Print["Quartic raw source decomposition from frozen W3"];
Print[""];
Print["Frozen cubic convention check H3 - {H2,W3} == 0: ",
  TrueQ[Simplify[Expand[H3 - bracketRaw[H2, W3]]] === 0]
];
Print[""];
Print["Explicit quartic source pieces:"];
Print["  S4H4 = ", S4H4];
Print[""];
Print["  S4W3Induced1Raw = {H3,W3} = ", S4W3Induced1Raw];
Print[""];
Print["  S4W3Induced2Raw = 1/2 {{H2,W3},W3} = ", S4W3Induced2Raw];
Print[""];
Print["  S4W3Induced = -{H3,W3} + 1/2 {{H2,W3},W3} = ", S4W3Induced];
Print[""];
Print["  S4Raw = S4H4 + S4W3Induced = ", S4Raw];
Print[""];
Print["Equivalence to current shorthand after using frozen W3 relation:"];
Print["  current shorthand H4Direct - 1/2 {H3,W3} = ", currentShorthand];
Print["  raw-minus-shorthand = ", equivalenceResidual];
Print["  check via H3 - {H2,W3}: ", Simplify[equivalenceViaFrozenW3]];
