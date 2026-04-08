(* Creamer 1989 tutorial, section 4.
   Goal:
     Explain why the 1D deep-water case becomes geometrically transparent.
*)

ClearAll["Global`*"];

$Assumptions = {
  k > 0,
  Element[{k}, Reals]
};

Print["Section 4: one-dimensional specialization and Hilbert geometry"];
Print[""];
Print["Equation (4.1): in one dimension the kernel D vanishes and B simplifies strongly."];
Print["This is the first sign that the 1D deep-water case is special."];
Print[""];
Print["Why does section 4 simplify so much?"];
Print["  Because in 1D the kernel structure collapses enough that the Lie-transform"];
Print["  equations can be rewritten using a Hilbert transform."];
Print[""];
Print["Equation (4.2): define the Hilbert transform"];
Print["  ZTilde(x,lambda) = (1/pi) P Integral dy Z(y,lambda)/(x-y)"];
Print["and similarly for PhiTilde."];
Print[""];
Print["Equation (4.3): in Fourier space"];
Print["  ZTilde_k = - i sgn(k) Z_k"];
Print["So for positive k, the Hilbert transform turns cosine into sine."];
Print[""];
Print["Paper interpretation right after (4.3):"];
Print["  in the linear approximation, ZTilde is the horizontal displacement of"];
Print["  fluid elements at the surface."];
Print[""];
Print["Why is that important?"];
Print["  Because the canonical transform becomes geometrically visible as a horizontal"];
Print["  rearrangement of the surface profile."];
Print[""];
Print["Equations (4.7a)-(4.7b): the 1D Lie-transform equations become"];
Print["  partial_lambda ZTilde   = -(1/2) partial_x (ZTilde^2)"];
Print["  partial_lambda PhiTilde = - ZTilde partial_x PhiTilde"];
Print[""];
Print["Equations (4.8)-(4.10): these are solved by characteristics."];
Print["The shifted coordinate chi is defined implicitly by"];
Print["  x = chi + lambda ZTilde_0(chi)"];
Print["or, when integrating backward from lambda = 1,"];
Print["  x = chi - (1-lambda) ZTilde_0(chi)."];
Print[""];
Print["Main takeaway:"];
Print["  The later geometric language of horizontal remapping comes from this"];
Print["  characteristic coordinate chi, not from an extra approximation added by hand."];
