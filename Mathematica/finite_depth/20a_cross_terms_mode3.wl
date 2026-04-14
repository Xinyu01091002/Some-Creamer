(* ::Package:: *)
(*  Cross terms at output mode 3, modes {-3,...,3}.
    All cross-term coefficients evaluated numerically at three kh values.
    T[n] = n*tanh(n*kh) is substituted as a pure number.
*)

ClearAll["Global`*"];

modes = Flatten[{Range[-3,-1], Range[1,3]}];

thetaOf[0] := 0;
thetaOf[n_Integer] /; n != 0 := T[Abs[n]];
omegaOf[n_Integer] /; n != 0 := Sqrt[T[Abs[n]]];
ksqOf[n_Integer] := n^2;

denBase[t1_,t2_,t3_] := t1(t2+t3-t1)+t2(t1+t3-t2)+t3(t1+t2-t3);

dFinite[t1_,t2_,t3_,q1_,q2_,q3_] := Module[{den,num},
  den = denBase[t1,t2,t3];
  If[t1===0||t2===0||t3===0||TrueQ[den==0],0,
    num = q1(t1^2-(t2-t3)^2)+q2(t2^2-(t1-t3)^2)+q3(t3^2-(t1-t2)^2)
        - 2 t1 t2 t3(t1+t2+t3);
    num/(12 den)]];

bFinite[t1_,t2_,t3_,q1_,q2_,q3_] := Module[{den,num},
  den = denBase[t1,t2,t3];
  If[TrueQ[den==0],0,
    num = t3(t3(t1+t2)-t1^2-t2^2)+t3(q1+q2-2q3)+(t1-t2)(q1-q2);
    num/(2 den)]];

DkerVal[a_,b_,c_] := DkerVal[a,b,c] =
  dFinite[thetaOf[a],thetaOf[b],thetaOf[c],ksqOf[a],ksqOf[b],ksqOf[c]];
BphiKerVal[a_,b_,c_] := BphiKerVal[a,b,c] =
  bFinite[thetaOf[a],thetaOf[b],thetaOf[c],ksqOf[a],ksqOf[b],ksqOf[c]];
H3kerVal[a_,b_,c_] := 1/4(ksqOf[b]+ksqOf[c]-ksqOf[a]-2 thetaOf[b] thetaOf[c]);

bracketRaw[F_,G_] := Sum[
  D[F,z[k]] D[G,p[-k]] - D[F,p[k]] D[G,z[-k]], {k,modes}];

H2 = 1/2 Sum[z[k] z[-k]+thetaOf[k] p[k] p[-k],{k,modes}];
H3 = Sum[If[MemberQ[modes,-(a+b)]&&a+b!=0,
    z[-(a+b)] p[a] p[b] H3kerVal[-(a+b),a,b],0],{a,modes},{b,modes}];
W3 = Sum[If[a+b+c==0,
    DkerVal[a,b,c] p[a] p[b] p[c]+BphiKerVal[a,b,c] z[a] z[b] p[c],0],
  {a,modes},{b,modes},{c,modes}];

N2PhiFn[k_] := N2PhiFn[k] = Sum[
  If[MemberQ[modes,d]&&MemberQ[modes,k-b-d],
    z[b] z[k-b-d] p[d](thetaOf[k] thetaOf[k-b] thetaOf[d]
      -1/2 thetaOf[k] ksqOf[d]-1/2 ksqOf[k] thetaOf[d]),0],
  {b,modes},{d,modes}];
S4H4 = 1/2 Sum[If[MemberQ[modes,-k],p[-k] N2PhiFn[k],0],{k,modes}];
S4W3 = Expand[-bracketRaw[H3,W3]+1/2 bracketRaw[bracketRaw[H2,W3],W3]];
S4Raw = Expand[S4H4+S4W3];

toN = Flatten[Table[{z[k]->(aa[k]+bb[k])/2,
  p[k]->(aa[k]-bb[k])/(2 I omegaOf[k])},{k,modes}]];
fromN = Flatten[Table[{aa[k]->z[k]+I omegaOf[k] p[k],
  bb[k]->z[k]-I omegaOf[k] p[k]},{k,modes}]];
normalVars = Flatten[Table[{aa[k],bb[k]},{k,modes}]];

nwt[mon_] := Module[{pa,pb},
  pa=Total[Table[Exponent[mon,aa[k]] omegaOf[k],{k,modes}]];
  pb=Total[Table[Exponent[mon,bb[k]] omegaOf[k],{k,modes}]];
  pa-pb];

splitNonres[poly_] := Module[{terms,mon,wt,nonres={}},
  terms=If[Head[Expand[poly]]===Plus,List@@Expand[poly],{Expand[poly]}];
  Do[mon=Times@@Table[var^Exponent[t,var],{var,normalVars}];
     wt=nwt[mon];
     If[!TrueQ[Simplify[wt==0]],AppendTo[nonres,t]],{t,terms}];
  Total[nonres]];

solveFn[poly_] := Module[{terms,mon,c,wt,solved},
  terms=If[Head[Expand[poly]]===Plus,List@@Expand[poly],{Expand[poly]}];
  solved=Table[
    mon=Times@@Table[var^Exponent[t,var],{var,normalVars}];
    c=Together[t/mon]; wt=nwt[mon];
    If[TrueQ[Simplify[wt==0]],0,-c mon/(I wt)],{t,terms}];
  Total[solved]];

W4H4 = Expand[solveFn[splitNonres[Expand[S4H4/.toN]]]/. fromN];
W4W3 = Expand[solveFn[splitNonres[Expand[S4W3/.toN]]]/. fromN];
W4Total = Expand[W4H4+W4W3];

Print["W3 and W4Total built. Building coordinate maps for mode 3 ..."];

dW3dp3   = D[W3, p[-3]];
d2W3dp3  = Sum[D[dW3dp3,z[k]] D[W3,p[-k]] - D[dW3dp3,p[k]] D[W3,z[-k]],{k,modes}];
dW4dp3   = D[W4Total, p[-3]];

mH3   = Expand[dW3dp3 + 1/2 d2W3dp3];
mFull = Expand[mH3 + dW4dp3];

Print["Maps built. Extracting cross-term expressions ..."];

(* wave state substitution: each mode gets its own symbol *)
(* z[k] = (A[|k|]/2), p[k] = sign * i * A[|k|] / (2 sqrt(T[|k|])) for k>0 *)
waveSub = Flatten[Table[{
  z[k]  -> A[Abs[k]]/2,
  z[-k] -> A[k]/2,
  p[k]  -> -I A[k]/(2 Sqrt[T[k]]),
  p[-k] ->  I A[k]/(2 Sqrt[T[k]])
},{k,1,3}]];

(* Substitute and expand in A[1], A[2], A[3] *)
mH3sub   = Expand[mH3   /. waveSub];
mFullsub = Expand[mFull /. waveSub];

(* Extract cubic coefficients in amplitude variables *)
getC[expr_,n1_,n2_,n3_] := 2 Coefficient[Coefficient[Coefficient[expr,A[1],n1],A[2],n2],A[3],n3];

(* Single-freq *)
c111H3   = getC[mH3sub,   3, 0, 0];
c111Full = getC[mFullsub, 3, 0, 0];

(* Bichromatic: need A[1]^1 A[2]^2 at mode 3 = (-1,2,2) triad *)
c122H3   = getC[mH3sub,   1, 2, 0];
c122Full = getC[mFullsub, 1, 2, 0];

(* Trichromatic: (-1,1,3)->3, (-2,2,3)->3 *)
c113H3   = getC[mH3sub,   2, 0, 1];
c113Full = getC[mFullsub, 2, 0, 1];

c223H3   = getC[mH3sub,   0, 2, 1];
c223Full = getC[mFullsub, 0, 2, 1];

Print["Coefficients extracted. Evaluating numerically ..."];
Print[""];

(* T[n] = n*tanh(n*kh) for any n, including intermediate wavenumbers from H4 *)
khList = {1/2, 1, 2, Infinity};
khLabels = {"kh=0.5", "kh=1", "kh=2", "DW"};

(* Cover T[1..9] to handle all intermediate modes in H4 kernel *)
tValsAt[kh_] := If[kh===Infinity,
  Table[T[n]->n, {n,1,9}],
  Table[T[n]->N[n Tanh[n kh], 20], {n,1,9}]];

evalAt[expr_, kh_] := N[expr /. tValsAt[kh], 10];

Print["=== A. Single-frequency sanity: (1,1,1)->3, eps1^3 ==="];
Print["  Verify H3-only -> H3-target, H3+H4 -> Stokes C33"];
Print["  Stokes C33 at kh=1: ", N[(27-9sigma^2+9sigma^4-3sigma^6)/(64sigma^6) /. sigma->Tanh[1], 10]];
Do[
  tv = tValsAt[kh];
  h3val   = N[c111H3   /. tv, 10];
  fullval = N[c111Full /. tv, 10];
  Print[StringForm["  `1`: H3-only = `2`, H3+H4 = `3`",
    khLabels[[idx]], Re[h3val], Re[fullval]]],
  {idx,1,4,1},{kh,{1/2,1,2,Infinity}}[[{idx}]]];
Print[""];

Print["=== B. Bichromatic cross: (-1,2,2)->3, eps1*eps2^2 ==="];
Print["  Label: H3-only, H3+H4, Delta, (H3+H4)/(H3-only)"];
Do[
  kh = {1/2,1,2,Infinity}[[idx]];
  tv = tValsAt[kh];
  h3val   = Re[N[c122H3   /. tv, 10]];
  fullval = Re[N[c122Full /. tv, 10]];
  delta   = fullval - h3val;
  ratio   = If[Abs[h3val]>10^-14, fullval/h3val, "---"];
  Print[StringForm["  `1`: H3=`2`, H3+H4=`3`, D=`4`, ratio=`5`",
    khLabels[[idx]], h3val, fullval, delta, ratio]],
  {idx,1,4}];
Print[""];

Print["=== C. Trichromatic cross: (-1,1,3)->3, eps1^2*eps3 ==="];
Do[
  kh = {1/2,1,2,Infinity}[[idx]];
  tv = tValsAt[kh];
  h3val   = Re[N[c113H3   /. tv, 10]];
  fullval = Re[N[c113Full /. tv, 10]];
  delta   = fullval - h3val;
  ratio   = If[Abs[h3val]>10^-14, fullval/h3val, "---"];
  Print[StringForm["  `1`: H3=`2`, H3+H4=`3`, D=`4`, ratio=`5`",
    khLabels[[idx]], h3val, fullval, delta, ratio]],
  {idx,1,4}];
Print[""];

Print["=== D. Trichromatic cross: (-2,2,3)->3, eps2^2*eps3 ==="];
Do[
  kh = {1/2,1,2,Infinity}[[idx]];
  tv = tValsAt[kh];
  h3val   = Re[N[c223H3   /. tv, 10]];
  fullval = Re[N[c223Full /. tv, 10]];
  delta   = fullval - h3val;
  ratio   = If[Abs[h3val]>10^-14, fullval/h3val, "---"];
  Print[StringForm["  `1`: H3=`2`, H3+H4=`3`, D=`4`, ratio=`5`",
    khLabels[[idx]], h3val, fullval, delta, ratio]],
  {idx,1,4}];
Print[""];

Print["=== Notes ==="];
Print["Coefficient C is defined so eta_nl[m] ~ (C/2) A1^n1 A2^n2 A3^n3 cos(m x)."];
Print["Factor of 2 already applied above. Imaginary parts should vanish (real waves)."];
Print["DW = deep water: T[1]=1, T[2]=2, T[3]=3."];
