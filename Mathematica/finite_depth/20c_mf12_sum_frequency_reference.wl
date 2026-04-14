(* ::Package:: *)
(*  MF12 sum-frequency reference coefficients, same getC convention as 20b.
    1D unidirectional, g=1, k1=1, k2=2, k3=3.
    kh sweeps {0.5, 1, 2}.

    Convention:
      eta_33^{nd} = A^3 * G3 / (2*(k1*h)^2) * cos(3x)  =>  getC_MF12 = G / (2*kh^2)
*)

ClearAll["Global`*"];

g = 1;
k1 = 1;  k2 = 2;  k3 = 3;

(* ------------------------------------------------------------------ *)
(*  Dispersion                                                         *)
(* ------------------------------------------------------------------ *)
omegaFn[kap_, h_] := If[h === Infinity, Sqrt[kap], Sqrt[kap Tanh[kap h]]];

(* ------------------------------------------------------------------ *)
(*  Second-order kernels (all args renamed to avoid pattern clash)     *)
(* ------------------------------------------------------------------ *)
Lambda2fn[won_, kn_, wom_, km_, womNpm_, aNpm_, gNpm_, bNpm_] :=
  Module[{knkm = kn km},
    1/(2 won wom bNpm) (
      g aNpm (won (km^2 + knkm) + wom (kn^2 + knkm))
      + gNpm (g^2 knkm + won^2 wom^2 - won wom womNpm^2))];

Gamma2fn[won_, kn_, wom_, km_, womNpm_, bNpm_] :=
  Module[{knkm = kn km},
    1/(2 won wom bNpm) (
      won wom womNpm (womNpm^2 - won wom)
      - g^2 won (km^2 + 2 knkm)
      - g^2 wom (kn^2 + 2 knkm))];

(* Note: Lambda2 and Gamma2 above absorbed one factor of h each since
   the original eq has h/(2 won wom beta).  But h appears in h-form in
   the original code.  Let me restore explicitly. *)

Lambda2fn[won_, kn_, wom_, km_, womNpm_, aNpm_, gNpm_, bNpm_, h_] :=
  Module[{knkm = kn km},
    h/(2 won wom bNpm) (
      g aNpm (won (km^2 + knkm) + wom (kn^2 + knkm))
      + gNpm (g^2 knkm + won^2 wom^2 - won wom womNpm^2))];

Gamma2fn[won_, kn_, wom_, km_, womNpm_, bNpm_, h_] :=
  Module[{knkm = kn km},
    h/(2 won wom bNpm) (
      won wom womNpm (womNpm^2 - won wom)
      - g^2 won (km^2 + 2 knkm)
      - g^2 wom (kn^2 + 2 knkm))];

(* Compute G_npm and F_npm for a sum pair kn+km (both positive, 1D) *)
npmPairFn[kn_, km_, h_] := Module[
  {won, wom, kNpm, wNpm, aNpm, gNpm, bNpm, GG, FF},
  won  = omegaFn[kn, h];
  wom  = omegaFn[km, h];
  kNpm = kn + km;
  wNpm = won + wom;
  aNpm = wNpm Cosh[kNpm h];
  gNpm = kNpm Sinh[kNpm h];
  bNpm = wNpm^2 Cosh[kNpm h] - g kNpm Sinh[kNpm h];
  GG   = Lambda2fn[won, kn, wom, km, wNpm, aNpm, gNpm, bNpm, h];
  FF   = Gamma2fn[won, kn, wom, km, wNpm, bNpm, h];
  {GG, FF, kNpm, gNpm, wNpm}];

(* ------------------------------------------------------------------ *)
(*  G_3: single-mode third-order surface kernel, Eq. 3.64             *)
(* ------------------------------------------------------------------ *)
G3fn[kap_, h_] :=
  (3 h^2 kap^2 / (128 Sinh[kap h]^6)) *
  (14 + 15 Cosh[2 kap h] + 6 Cosh[4 kap h] + Cosh[6 kap h]);

(* ------------------------------------------------------------------ *)
(*  Lambda3: Eq. 3.53                                                  *)
(*  Arguments: three modes (won,kn), (wom,km), (wop,kp)              *)
(*             three second-order pairs: nm, np, mp                   *)
(*             output mode: wout, aout, gout, bout                    *)
(* ------------------------------------------------------------------ *)
Lambda3fn[won_,kn_, wom_,km_, wop_,kp_,
          kNm_,gNm_,GNm_,FNm_,
          kNp_,gNp_,GNp_,FNp_,
          kMp_,gMp_,GMp_,FMp_,
          wout_,aout_,gout_,bout_,h_] :=
  Module[{knkm = kn km, knkp = kn kp, kmkp = km kp},
  h^2/(4 bout) (
    aout (won(knkm+knkp+kn^2) + wom(knkm+kmkp+km^2) + wop(knkp+kmkp+kp^2))
    + gout (
        g/won (wom knkm + wop knkp - wout kn^2)
      + g/wom (won knkm + wop kmkp - wout km^2)
      + g/wop (won knkp + wom kmkp - wout kp^2)))
  - h FNm/(2 bout) (
      aout Cosh[h kNm] (knkp+kmkp+kNm^2)
    + gout (g/wop (knkp+kmkp) Cosh[h kNm] - gNm wout))
  - h FNp/(2 bout) (
      aout Cosh[h kNp] (knkm+kmkp+kNp^2)
    + gout (g/wom (knkm+kmkp) Cosh[h kNp] - gNp wout))
  - h FMp/(2 bout) (
      aout Cosh[h kMp] (knkm+knkp+kMp^2)
    + gout (g/won (knkm+knkp) Cosh[h kMp] - gMp wout))
  + h GNm/(2 bout) (aout g/wop (knkp+kmkp+kp^2) - gout wop^2)
  + h GNp/(2 bout) (aout g/wom (knkm+kmkp+km^2) - gout wom^2)
  + h GMp/(2 bout) (aout g/won (knkm+knkp+kn^2) - gout won^2)];

(* Compute G for triad (kn,km,kp) -> kout = kn+km+kp, all positive 1D *)
triG[kn_, km_, kp_, h_] := Module[
  {won, wom, wop,
   pnm, pnp, pmp, GNm, FNm, kNm, gNm, GNp, FNp, kNp, gNp, GMp, FMp, kMp, gMp,
   kout, wout, aout, gout, bout},
  won = omegaFn[kn, h];
  wom = omegaFn[km, h];
  wop = omegaFn[kp, h];
  pnm = npmPairFn[kn, km, h];  {GNm, FNm, kNm, gNm, _} = pnm;
  pnp = npmPairFn[kn, kp, h];  {GNp, FNp, kNp, gNp, _} = pnp;
  pmp = npmPairFn[km, kp, h];  {GMp, FMp, kMp, gMp, _} = pmp;
  kout = kn + km + kp;
  wout = won + wom + wop;
  aout = wout Cosh[kout h];
  gout = kout Sinh[kout h];
  bout = wout^2 Cosh[kout h] - g kout Sinh[kout h];
  Lambda3fn[won,kn, wom,km, wop,kp,
            kNm,gNm,GNm,FNm, kNp,gNp,GNp,FNp, kMp,gMp,GMp,FMp,
            wout,aout,gout,bout,h]];

(* ------------------------------------------------------------------ *)
(*  Convert G to getC convention: getC_MF12 = G / (2 * (k_ref*h)^2)  *)
(*  k_ref = k1 = 1, so getC = G / (2*kh^2)                           *)
(* ------------------------------------------------------------------ *)
toGetC[GG_, kh_] := GG / (2 kh^2);

(* ------------------------------------------------------------------ *)
(*  Creamer H3+H4 results from 20b (hardcoded for comparison table)   *)
(* ------------------------------------------------------------------ *)
(* rows: kh=0.5, kh=1, kh=2 *)
crH3H4 = {
  (* (1,1,1)->3 *) {40.846, 1.9395, 0.4672},
  (* (1,1,2)->4 *) {67.891, 5.8878, 2.2904},
  (* (1,2,2)->5 *) {38.726, 5.8052, 3.3382},
  (* (2,2,2)->6 *) {7.758,  1.8688, 1.5061},
  (* (1,1,3)->5 *) {60.988, 8.0280, 3.5453},
  (* (1,2,3)->6 *) {70.433, 15.224, 9.5739}
};
triadLabels = {"(1,1,1)->3", "(1,1,2)->4", "(1,2,2)->5",
               "(2,2,2)->6", "(1,1,3)->5", "(1,2,3)->6"};

(* ------------------------------------------------------------------ *)
(*  Main loop                                                          *)
(* ------------------------------------------------------------------ *)
khValues = {1/2, 1, 2};
khLabels = {"kh=0.5", "kh=1", "kh=2"};

(* Compute MF12 getC for all triads at each kh *)
(*  Normalization notes
    --------------------
    MF12 ThetaA always divides by h^2 (physical depth, same for all modes).
    So the reference depth is always k1*h = kh, not mode-specific k_n*h.

    For self-interaction / repeated-mode triads (G_3, G_np2m, G_2npm):
      eta = A_triad * G * cos(m x)  with  A_triad = a1*...*an/(2*h^2)
      => getC_MF12 = G / (2*kh^2)  = toGetC[G, kh]

    For distinct-mode triplets (G_npmpp, n,m,p all different):
      MF12 code uses  Z = 2 * A_npmpp  (permutation symmetry factor)
      => eta = a1*a2*a3/h^2 * G_npmpp * cos(m x)
      => getC_MF12 = G_npmpp / kh^2  = 2 * toGetC[G_npmpp, kh]
*)
mf12Table = Table[
  kh = N[khValues[[j]], 20];
  h  = kh;
  {
    toGetC[G3fn[k1, h], kh],          (* (1,1,1)->3: self-interact at k1 *)
    toGetC[triG[k1,k1,k2, h], kh],   (* (1,1,2)->4: G_2npm, 2*k1+k2 *)
    toGetC[triG[k1,k2,k2, h], kh],   (* (1,2,2)->5: G_np2m, k1+2*k2 *)
    toGetC[G3fn[k2, h], kh],          (* (2,2,2)->6: self-interact at k2, ref depth = k1*h *)
    toGetC[triG[k1,k1,k3, h], kh],   (* (1,1,3)->5: G_2npm with k3 *)
    2 toGetC[triG[k1,k2,k3, h], kh]  (* (1,2,3)->6: G_npmpp, distinct modes -> factor 2 *)
  },
  {j, 1, Length[khValues]}];

Print["=== MF12 sum-frequency getC coefficients ==="];
Print["getC_MF12 = G/(2*kh^2),  same convention as 20b Creamer output"];
Print[""];
Do[
  kh = N[khValues[[j]], 6];
  Print["--- " <> khLabels[[j]] <> " ---"];
  Do[
    mfv = N[mf12Table[[j, i]], 8];
    Print["  " <> triadLabels[[i]] <> "  getC_MF12 = " <>
          ToString[NumberForm[mfv, {10,5}]]],
    {i, 1, 6}];
  Print[""],
  {j, 1, 3}];

Print["=== Comparison table: H3+H4 / MF12 ratio ==="];
Print["ratio > 1 : Creamer overshoots MF12"];
Print["ratio < 1 : Creamer undershoots MF12 (gap remains)"];
Print["ratio = 1 : exact agreement"];
Print[""];
Print["Triad        | kh=0.5  H3+H4  MF12  ratio | kh=1    H3+H4  MF12  ratio | kh=2   H3+H4  MF12  ratio"];
Do[
  row = triadLabels[[i]] <> " |";
  Do[
    mfv = N[mf12Table[[j, i]], 6];
    crv = crH3H4[[i, j]];
    rat = N[crv/mfv, 5];
    row = row <> " " <> ToString[NumberForm[mfv,{7,3}]] <>
                 " " <> ToString[NumberForm[crv,{7,3}]] <>
                 " " <> ToString[NumberForm[rat,{5,3}]] <> " |",
    {j, 1, 3}];
  Print[row],
  {i, 1, 6}];

Print[""];
Print["Deep water: all H3+H4 = MF12 = Stokes (ratio = 1 exactly, see 20b)."];
