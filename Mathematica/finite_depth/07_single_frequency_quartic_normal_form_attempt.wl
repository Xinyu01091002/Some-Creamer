(* ::Package:: *)
(*
  Attempt a direct single-frequency quartic normal-form calculation.

  This script is an executable algebra scaffold for the hard step:

      K4 = H4 + 1/2 {H3, W3}

  in a finite set of 1-D Fourier modes.  It uses the canonical bracket

      {F,G} = Sum_k dF/dz_k dG/dp_-k - dF/dp_k dG/dz_-k

  and the finite-depth kernels already used by the H3-only lambda-flow.

  The current purpose is to compute and expose the one-frequency quartic
  Hamiltonian data.  A complete W4(k1,k2,k3,k4) solve will use this scaffold,
  but this script deliberately prints the obstruction/target rather than
  silently assuming a W4 coefficient.
*)

ClearAll["Global`*"];

$Assumptions = Element[sigma, Reals] && 0 < sigma < 1;

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

Bker[a_, b_, c_] := Bker[a, b, c] =
  bFinite[theta[a], theta[b], theta[c], ksq[a], ksq[b], ksq[c]];

H3ker[a_, b_, c_] :=
  1/4 (ksq[b] + ksq[c] - ksq[a] - 2 theta[b] theta[c]);

(* Ordered H4 kernel from 04:
   outer phi mode a, eta modes b,c, inner phi mode d. *)
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

rulesSigma = {
  T[1] -> sigma,
  T[2] -> 2 Tanh[2 ArcTanh[sigma]],
  T[3] -> 3 Tanh[3 ArcTanh[sigma]],
  T[4] -> 4 Tanh[4 ArcTanh[sigma]],
  T[5] -> 5 Tanh[5 ArcTanh[sigma]],
  T[6] -> 6 Tanh[6 ArcTanh[sigma]]
};

monomial = z[1]^3 p[-3];

h4CoeffRaw = Coefficient[H4, monomial];
bracketCoeffRaw = Coefficient[bracketRaw[H3, W3], monomial];
k4CoeffRaw = Together[h4CoeffRaw + 1/2 bracketCoeffRaw];

h4CoeffSigma = FullSimplify[h4CoeffRaw /. rulesSigma];
bracketCoeffSigma = FullSimplify[bracketCoeffRaw /. rulesSigma];
k4CoeffSigma = FullSimplify[k4CoeffRaw /. rulesSigma];

requiredW4SurfaceCoeff =
  (3 (-1 + sigma^2)^2 (27 + 60 sigma^2 + 50 sigma^4 + 4 sigma^6 + 3 sigma^8))/
    (16 sigma^6 (9 + 14 sigma^2 + 9 sigma^4));

Print["Single-frequency quartic normal-form scaffold"];
Print[""];
Print["Computed K4 = H4 + 1/2 {H3,W3} in modes +/-1,+/-2,+/-3."];
Print["Coefficient contributions for monomial z[1]^3 p[-3]:"];
Print["  H4 contribution:       ", h4CoeffSigma];
Print["  {H3,W3} contribution:  ", bracketCoeffSigma];
Print["Coefficient of monomial z[1]^3 p[-3] in K4:"];
Print["  raw T[n] form: ", k4CoeffRaw];
Print["  sigma form:    ", k4CoeffSigma];
Print[""];
Print["Previously calibrated minimal one-mode W4 surface coefficient w:"];
Print["  ", requiredW4SurfaceCoeff];
Print[""];
Print["Deep-water limit of K4 monomial coefficient:"];
Print["  ", FullSimplify[k4CoeffSigma /. sigma -> 1]];
Print["Deep-water limit of calibrated w:"];
Print["  ", FullSimplify[requiredW4SurfaceCoeff /. sigma -> 1]];
Print[""];
Print["Status:"];
Print["  This script constructs K4 explicitly in a truncated mode algebra."];
Print["  The remaining step is to apply the homological inverse of ad_H2 in"];
Print["  diagonal normal variables and map the resulting W4 back to zeta."];
Print["  Do not equate the K4 monomial coefficient directly with the surface"];
Print["  W4 coefficient; they live in different parts of the normal-form chain."];
