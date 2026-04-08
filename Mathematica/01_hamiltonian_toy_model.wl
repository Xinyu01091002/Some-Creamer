(* Creamer 1989 tutorial, section 2.
   Goal:
     Explain why the full Hamiltonian H(zeta,phi_s) can be written as
       H = H2 + H3 + H4 + ...
     and what H2 and H3 mean in the paper.

   Important rule:
     No toy canonical pair is introduced here.
     This script stays at the paper level.
*)

ClearAll["Global`*"];

$Assumptions = {
  g > 0,
  Element[{g}, Reals]
};

Print["Section 2: from H to H2, H3, H4, ..."];
Print[""];
Print["Step 1. What are the basic variables?"];
Print["  zeta(x,t)   = free-surface elevation"];
Print["  phi_s(x,t)  = velocity potential evaluated at the free surface"];
Print[""];
Print["Step 2. Why do the operators D_x and D_z depend on zeta?"];
Print["  Because they are defined by first solving Laplace's equation in the fluid"];
Print["  interior with boundary data prescribed on the ACTUAL free surface"];
Print["  z = zeta(x,t)."];
Print["  If zeta changes, the boundary where the harmonic extension is matched changes,"];
Print["  so the surface operators change as well."];
Print[""];
Print["Step 3. Full Hamiltonian, equation (2.9):"];
Print["  H(zeta,phi_s) = (1/2) Integral d^2x"];
Print["                 [ phi_s ( D_z - (partial_x zeta) . D_x ) phi_s + g zeta^2 ]"];
Print[""];
Print["Step 4. Why can H be written as H2 + H3 + H4 + ... ?"];
Print["  Because D_x and D_z are zeta-dependent operators."];
Print["  In weakly nonlinear theory one expands about the flat reference surface z = 0."];
Print["  Then the operators themselves admit an expansion of the form"];
Print["    D = D^(0) + D^(1)[zeta] + D^(2)[zeta] + ..."];
Print[""];
Print["  Substituting that operator expansion into H produces terms of different"];
Print["  total degree in the canonical variables zeta and phi_s."];
Print[""];
Print["Step 5. Meaning of the notation"];
Print["  H2 = all terms of total degree 2 in (zeta,phi_s)"];
Print["  H3 = all terms of total degree 3"];
Print["  H4 = all terms of total degree 4"];
Print["  and so on."];
Print[""];
Print["Why is this degree counting reasonable?"];
Print["  Because weakly nonlinear theory treats zeta and phi_s as small-amplitude"];
Print["  variables of comparable perturbative order."];
Print[""];
Print["Step 6. Quadratic Hamiltonian, equation (2.10):"];
Print["  H2 = (1/2) Integral d^2x [ phi_s theta phi_s + g zeta^2 ]"];
Print["  with theta = sqrt( - (partial_x)^2 )"];
Print["  and in Fourier space theta(k) = |k|."];
Print[""];
Print["Why is H2 the linear Hamiltonian?"];
Print["  Because it is exactly the flat-surface part of the energy, with no nonlinear"];
Print["  correction from the zeta-dependence of the surface operators."];
Print[""];
Print["Step 7. Cubic Hamiltonian, equation (2.13):"];
Print["  H3 is the first nonlinear correction."];
Print["  In Fourier space it has the structure"];
Print["    zeta_1 phi_2 phi_3 times a k-dependent kernel,"];
Print["  together with a delta function enforcing k1 + k2 + k3 = 0."];
Print[""];
Print["Important clarification about the delta function:"];
Print["  The delta^2 in equation (2.13) enforces wavevector conservation only."];
Print["  It does NOT say that the interacting triad is resonant."];
Print["  Resonance would require, in addition, a frequency condition of the form"];
Print["    +/- omega1 +/- omega2 +/- omega3 = 0."];
Print["  For deep-water gravity waves the wavevector condition can hold while the"];
Print["  frequency condition still fails, so H3 can exist and yet remain non-resonant."];
Print[""];
Print["Why does H3 look like zeta * phi * phi ?"];
Print["  Because the first correction to the operators contributes one power of zeta,"];
Print["  while the Hamiltonian remains quadratic in phi_s outside those operators."];
Print["  So the first nonlinear correction naturally has one zeta factor and two phi_s factors."];
Print[""];
Print["Main takeaway:"];
Print["  H2 governs linear deep-water gravity waves."];
Print["  H3 is the leading nonlinear term that section 3 tries to remove by a"];
Print["  canonical transformation."];
