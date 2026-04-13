(* ::Package:: *)
(*
  Physical-scaling check for the direct H4 normal-form calculation.

  Script 11 works in the dimensionless convention k0 = g = 1 and compares the
  C33 coefficient.  This script keeps k0 and g symbolic/numeric in the same
  algebra used by the MATLAB 1D triad prototype, then evaluates the single-mode
  eta_hat(3k0) increment directly.  Its purpose is to locate the missing scaling
  between the Mathematica normal-form check and physical FFT coefficients.
*)

ClearAll["Global`*"];

modes = {-3, -2, -1, 1, 2, 3};

g0 = 9.81;
k0 = 0.0279;
kh0 = 2.0;
sigma0 = N[Tanh[kh0]];
eps0 = 0.05;
etaHat1 = eps0/(2 k0);

theta[0] := 0;
theta[n_Integer] /; n != 0 := Abs[n] k0 Tanh[Abs[n] kh0];
omega[n_Integer] /; n != 0 := Sqrt[g0 theta[n]];
gamma[n_Integer] /; n != 0 := Sqrt[theta[n]/g0];
ksq[n_Integer] := (n k0)^2;

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
    num/(12 g0 den)
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
      p[k] -> (aa[k] - bb[k])/(2 I gamma[k])
    },
    {k, modes}
  ]];

fromNormalRules = Flatten[Table[
    {
      aa[k] -> z[k] + I gamma[k] p[k],
      bb[k] -> z[k] - I gamma[k] p[k]
    },
    {k, modes}
  ]];

normalVars = Flatten[Table[{aa[k], bb[k]}, {k, modes}]];

normalWeight[mon_] := Module[{pa, pb},
  pa = Total[Table[Exponent[mon, aa[k]] omega[k], {k, modes}]];
  pb = Total[Table[Exponent[mon, bb[k]] omega[k], {k, modes}]];
  pa - pb
];

solveHomologicalNonres[poly_] := Module[{terms, mon, coeff, weight, solved},
  terms = If[Head[Expand[poly]] === Plus, List @@ Expand[poly], {Expand[poly]}];
  solved = Table[
    mon = Times @@ Table[var^Exponent[term, var], {var, normalVars}];
    coeff = Together[term/mon];
    weight = normalWeight[mon];
    If[TrueQ[FullSimplify[weight == 0]], 0, -coeff mon/(I weight)],
    {term, terms}
  ];
  Total[solved]
];

linearWaveRules = Flatten[Table[
    {
      z[k] -> Which[k == 1, etaHat1, k == -1, etaHat1, True, 0],
      p[k] -> Which[
        k == 1, -I Sqrt[g0/theta[1]] etaHat1,
        k == -1, I Sqrt[g0/theta[1]] etaHat1,
        True, 0
      ]
    },
    {k, modes}
  ]];

H4Normal = Expand[H4Direct /. toNormalRules];
BracketNormal = Expand[bracketRaw[H3, W3] /. toNormalRules];
K4 = Expand[H4Normal - 1/2 BracketNormal];
W4Normal = solveHomologicalNonres[K4];
W4Back = Expand[W4Normal /. fromNormalRules];

targetMonA = aa[-3] aa[1]^3;
targetMonB = bb[-3] aa[1]^3;
targetK4A = Coefficient[K4, targetMonA];
targetK4B = Coefficient[K4, targetMonB];
targetH4A = Coefficient[H4Normal, targetMonA];
targetH4B = Coefficient[H4Normal, targetMonB];
targetBracketA = Coefficient[BracketNormal, targetMonA];
targetBracketB = Coefficient[BracketNormal, targetMonB];
targetW4A = Coefficient[W4Normal, targetMonA];
targetW4B = Coefficient[W4Normal, targetMonB];

eta33H3Hat = Chop[N[(1/2 bracketRaw[bracketRaw[z[3], W3], W3]) /. linearWaveRules]];
eta33W4Hat = Chop[N[bracketRaw[z[3], W4Back] /. linearWaveRules]];

c33H3 = eta33H3Hat/(4 k0^2 etaHat1^3);
c33W4 = eta33W4Hat/(4 k0^2 etaHat1^3);
c33Stokes = (27 - 9 sigma0^2 + 9 sigma0^4 - 3 sigma0^6)/(64 sigma0^6);
c33H3Formula =
  (3*(1 + sigma0^2)*(54 + 39*sigma0^2 + 70*sigma0^4 - 28*sigma0^6 - 4*sigma0^8 - 3*sigma0^10))/
    (64*sigma0^6*(9 + 14*sigma0^2 + 9*sigma0^4));

Print["Physical scaling H4 check"];
Print["  g0 = ", N[g0, 12], ", k0 = ", N[k0, 12], ", kh0 = ", kh0,
  ", sigma = ", N[sigma0, 12], ", eps = ", N[eps0, 12]];
Print["  etaHat1 = ", N[etaHat1, 12]];
Print[""];
Print["  eta33 H3 hat =       ", N[eta33H3Hat, 16]];
Print["  eta33 W4 hat =       ", N[eta33W4Hat, 16]];
Print["  eta33 H3+W4 hat =    ", N[eta33H3Hat + eta33W4Hat, 16]];
Print["  expected Stokes hat = ", N[4 c33Stokes k0^2 etaHat1^3, 16]];
Print[""];
Print["  K4 coeff aa[-3] aa[1]^3 = ", N[targetK4A, 16]];
Print["  K4 coeff bb[-3] aa[1]^3 = ", N[targetK4B, 16]];
Print["  H4 coeff aa[-3] aa[1]^3 = ", N[targetH4A, 16]];
Print["  H4 coeff bb[-3] aa[1]^3 = ", N[targetH4B, 16]];
Print["  bracket coeff aa[-3] aa[1]^3 = ", N[targetBracketA, 16]];
Print["  bracket coeff bb[-3] aa[1]^3 = ", N[targetBracketB, 16]];
Print["  W4 coeff aa[-3] aa[1]^3 = ", N[targetW4A, 16]];
Print["  W4 coeff bb[-3] aa[1]^3 = ", N[targetW4B, 16]];
Print[""];
Print["  C33 H3 from bracket = ", N[c33H3, 16]];
Print["  C33 H3 formula =      ", N[c33H3Formula, 16]];
Print["  C33 W4 from bracket = ", N[c33W4, 16]];
Print["  C33 required delta =  ", N[c33Stokes - c33H3Formula, 16]];
Print["  C33 Stokes =          ", N[c33Stokes, 16]];
