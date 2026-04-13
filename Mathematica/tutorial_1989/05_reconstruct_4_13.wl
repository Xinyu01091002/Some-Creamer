(* Creamer 1989 tutorial, section 4 reconstruction.
   Goal:
     Present the paper-level reconstruction formulas faithfully and explain why
     they generate bound harmonics.
*)

ClearAll["Global`*"];

Print["Section 4 reconstruction: from linear variables back to the physical surface"];
Print[""];
Print["Equations (4.11)-(4.12): first write the linear variables in terms of the"];
Print["physical ones by using the inverse Hilbert transform and then shifting the"];
Print["integration variable to the characteristic coordinate chi."];
Print[""];
Print["Equations (4.13a)-(4.13b): physical variables in terms of the linear variables"];
Print["  zeta(x) = (1/pi) P Integral dy [1 - zetaTilde0'(y)] / [x - y + zetaTilde0(y)] zetaTilde0(y)"];
Print["  phi_s(x) is written analogously in terms of phiTilde_s0(y)."];
Print[""];
Print["This is the central reconstruction formula of stage 1."];
Print["It says: evolve the transformed variables linearly, then reconstruct the"];
Print["physical surface nonlinearly from them."];
Print[""];
Print["Why does this create bound harmonics?"];
Print["  Because the reconstruction is nonlinear in the parent field zetaTilde0."];
Print["  Nonlinearity enters through the map back to the physical surface, not through"];
Print["  introducing independent new free-wave modes."];
Print[""];
Print["Equation (4.14a): Fourier-space reconstruction for zeta_k"];
Print["  zeta_k = (1/|k|) Integral dy exp[-i k y] ( exp[i k zetaTilde0(y)] - 1 )"];
Print["This is one of the cleanest places to see why bound harmonics appear."];
Print["The exponential generates an infinite harmonic hierarchy from the parent field."];
Print[""];
Print["Equation (4.15): once reconstruction is known, spatial derivatives of the"];
Print["physical fields can also be written directly in terms of the linear variables."];
Print[""];
Print["Main takeaway:"];
Print["  The transformed variables evolve linearly, but the physical surface can still"];
Print["  be strongly non-sinusoidal because the reconstruction map is nonlinear."];
