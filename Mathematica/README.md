# Mathematica Workspace

This folder is reserved for Wolfram Mathematica and `wolframscript` derivations related to the Creamer transform project.

## Purpose

Use this directory for symbolic work that is awkward to do by hand, especially:

- series expansions
- Fourier-side algebra
- operator identities in 1-D
- canonical-transform bookkeeping
- checks of finite-depth asymptotics
- comparison of reconstructed wave shapes against benchmark formulas

## Current Status

- `wolframscript` is available on this machine and can be called from the terminal.
- Mathematica should be useful for expansion, simplification, and symbolic consistency checks.
- It may still be awkward for fully automated Lie-transform derivations unless the workflow is broken into smaller symbolic tasks.

## Layout

- `tutorial_1989/`
  - stage-1 scripts for understanding Creamer et al. (1989)
  - includes notation setup, toy Hamiltonian checks, Poisson/Lie tools, the 1-D Hilbert remap note, and the `(4.13)` reconstruction note
- `finite_depth/`
  - scripts for Wright-Creamer finite-depth questions
  - currently includes `01_single_frequency_limit.wl`, which derives the 1-D single-frequency small-steepness coefficients from the finite-depth canonical-pair lambda-flow
  - includes `02_hamiltonian_h2_h3_check.wl`, which checks the finite-depth `H2/H3` Hamiltonian pieces from the Dirichlet-Neumann expansion before attempting `H4/W4`
  - includes `03_creamer_1989_h2_h3_from_d_operators.wl`, which re-derives Creamer et al. (1989) `H2/H3` directly from the paper's `D_z - (partial_x zeta).D_x` notation
  - includes `04_hamiltonian_h4_from_d_operators.wl`, which derives the operator-space `N2` and `H4` source from the same D-operator expansion
  - includes `05_lie_h3_h4_single_frequency_check.wl`, which records the single-frequency `C[3,3]` correction a future `W4` calculation must supply
  - includes `06_single_frequency_w4_calibrated_stokes_check.wl`, which applies that calibrated one-mode correction and checks it exactly matches the finite-depth Stokes/VWA third-harmonic coefficient
  - includes `07_single_frequency_quartic_normal_form_attempt.wl`, which builds the truncated single-frequency `K4 = H4 + 1/2 {H3,W3}` algebra and exposes the target monomial coefficient before diagonal homological inversion
  - includes `08_single_frequency_w4_homological_inverse.wl`, which diagonalizes `H2`, solves the truncated homological inverse for `W4`, maps back to `(zeta, phi_s)`, and records the residual against the finite-depth Stokes/VWA `C[3,3]` target
  - includes `09_overall_h3_h4_absorption_single_frequency.wl`, which checks the full cubic-order coordinate map `zeta + {zeta,W3} + 1/2 {{zeta,W3},W3} + {zeta,W4}` and diagnoses whether the Hamiltonian-normal-form `W3` is aligned with the existing Creamer lambda-flow coordinates
  - includes `10_lambda_flow_w3_alignment_debug.wl`, which verifies directly that the `W3` bracket vector field reproduces the `01` lambda-flow H3-only `C[3,3]` result before adding quartic corrections
  - includes `11_direct_h4_convolution_overall_check.wl`, which rebuilds `H4` directly from the `N2` operator convolution, absorbs only the non-resonant quartic part, and checks the full single-frequency `C[3,3]` correction
  - includes `12_physical_scaling_h4_check.wl`, which converts the single-frequency quartic check to physical scaling and verifies the finite-depth `eta33` correction numerically against the supplied Stokes/VWA target
  - includes `13_h3_h4_overall_homological_check.wl`, which checks the overall finite-depth H3/H4 homological equations, separates resonant and non-resonant quartic parts, and compares single-frequency and toy-broadband observable residuals under the same `W3/W4` convention
  - includes `14_iterative_transform_search.wl`, which iterates through the most plausible non-resonant quartic sign conventions and adds checkpoints so the heavy exact residual search can be debugged stage by stage
  - includes `15_fast_numeric_transform_screen.wl`, which keeps the exact single-frequency check but uses sampled Hamiltonian residuals to quickly screen which quartic transform convention is worth a heavier exact proof

## Recommended Workflow

1. Start from a clearly stated symbolic target.
   - example: derive second-order bound harmonics from a horizontal remapping
   - example: compare deep-water and finite-depth kernels
2. Encode assumptions explicitly.
   - small steepness
   - real-valued fields
   - one-dimensional geometry unless otherwise stated
3. Use Mathematica for local algebraic tasks, not as a black-box theory engine.
4. Save each derivation as either:
   - a `.wl` script for reproducible terminal execution
   - a `.nb` notebook for exploratory work
5. Summarize successful derivations back into the repo notes.

## Suggested First Targets

- Reproduce the small-amplitude harmonic generation from a horizontal remapping.
- Build a symbolic comparison between deep-water and finite-depth linear operators.
- Test how far normal-form or Lie-transform manipulations can be pushed before manual intervention becomes necessary.

## Current Finite-Depth Script

Run:

```powershell
wolframscript -file Mathematica\finite_depth\01_single_frequency_limit.wl
```

The script uses `T[n] = n Tanh[n k0 h]` as the finite-depth symbol and keeps the default expansion to third order. Its deep-water check substitutes `T[n] -> n` and verifies that the leading single-frequency harmonic coefficients reduce to Taylor's deep-water Creamer/Stokes-leading values:

```text
1, 1/2, 3/8
```

Before working on quartic corrections, run:

```powershell
wolframscript -file Mathematica\finite_depth\02_hamiltonian_h2_h3_check.wl
```

For a more literal check against Creamer et al. (1989), equation `(2.9)`, run:

```powershell
wolframscript -file Mathematica\finite_depth\03_creamer_1989_h2_h3_from_d_operators.wl
```

This starts from the paper's `D_z - (partial_x zeta).D_x` operator, expands the interior harmonic potential around the flat surface, and verifies that the result gives Creamer et al. (1989), equations `(2.10)` and `(2.13)`.

The related `02` script starts from the finite-depth Dirichlet-Neumann form:

```text
G(eta) = Dz - Grad[eta].Dx
H = 1/2 Integral [g eta^2 + phi G(eta) phi] dx
G0(k) = |k| tanh(|k| h)
G1(eta) phi = -G0 eta G0 phi - Div[eta Grad[phi]]
```

It introduces a bookkeeping parameter `alpha`, expands `H(alpha eta, alpha phi)` by degree, identifies `H2` and `H3`, and then verifies that the resulting `H3` Fourier coefficient matches Wright-Creamer (1994), equation `(10)`, with the deep-water limit matching Creamer et al. (1989), equation `(2.13)`.

For the current quartic work, run:

```powershell
wolframscript -file Mathematica\finite_depth\04_hamiltonian_h4_from_d_operators.wl
wolframscript -file Mathematica\finite_depth\05_lie_h3_h4_single_frequency_check.wl
wolframscript -file Mathematica\finite_depth\06_single_frequency_w4_calibrated_stokes_check.wl
wolframscript -file Mathematica\finite_depth\07_single_frequency_quartic_normal_form_attempt.wl
wolframscript -file Mathematica\finite_depth\08_single_frequency_w4_homological_inverse.wl
wolframscript -file Mathematica\finite_depth\09_overall_h3_h4_absorption_single_frequency.wl
wolframscript -file Mathematica\finite_depth\10_lambda_flow_w3_alignment_debug.wl
wolframscript -file Mathematica\finite_depth\11_direct_h4_convolution_overall_check.wl
wolframscript -file Mathematica\finite_depth\12_physical_scaling_h4_check.wl
wolframscript -file Mathematica\finite_depth\13_h3_h4_overall_homological_check.wl
wolframscript -file Mathematica\finite_depth\14_iterative_transform_search.wl
wolframscript -file Mathematica\finite_depth\15_fast_numeric_transform_screen.wl
wolframscript -file Mathematica\finite_depth\16_frozen_w3_basis_check.wl
wolframscript -file Mathematica\finite_depth\17_quartic_raw_source_decomposition.wl
wolframscript -file Mathematica\finite_depth\18_split_nonres_quartic_solve.wl
wolframscript -file Mathematica\finite_depth\19_quartic_validation_from_frozen_w3.wl
```

The `04` script derives:

```text
N2 psi = G0 eta G0 eta G0 psi
         + 1/2 G0(eta^2 Delta psi)
         + 1/2 Delta(eta^2 G0 psi)
H4 = 1/2 Integral phi_s N2 phi_s dx
```

The `05` script is a diagnostic target, not a completed `W4` derivation. It records the finite-depth single-frequency `cos(3 theta)` correction required to make the current `H3`-only Creamer result match the supplied Stokes/VWA coefficient.

The `06` script applies this one-mode calibrated correction as an acceptance-test target. It also records the equivalent minimal one-mode generator coefficient for an ansatz of the form `W4_min = w z_1^3 phi_-3 + c.c.`. It verifies:

```text
Creamer_H3_only_C33 + calibrated_W4_delta_C33 - Stokes_VWA_C33 == 0
```

This is not yet a derivation of the full quartic generator `W4`; it is the single-frequency result that the later `K4 = H4 + 1/2 {H3,W3}` calculation must reproduce.

The `07` and `08` scripts are the first direct quartic normal-form attempt.  They verify that the current `W3` satisfies the cubic homological check on a sample triad, then compute the truncated `W4` contribution to the `z[1]^3 phi[-3]` surface coefficient.  With the current conventions this direct `W4` attempt does **not** reproduce the calibrated finite-depth Stokes/VWA correction; the printed residual is the next object to diagnose before claiming that absorbing `H4` is sufficient.

The `10` script corrects an earlier diagnostic worry: the `W3` bracket vector field does reproduce the existing `01` Creamer lambda-flow `C[3,3]` coefficient when the coordinate contribution `1/2 {{zeta,W3},W3}` is evaluated directly.  The `09` script now uses this overall coordinate-map check and explicitly splits `K4` into resonant and non-resonant monomials in normal variables.  Only the non-resonant quartic terms are passed to the `W4` homological inverse; resonant terms remain in the quartic normal form.  With the current `K4` and coordinate-map conventions, this non-resonant `W4` scan still does not exactly reproduce the finite-depth Stokes/VWA `C[3,3]` target, so the remaining issue is likely the precise quartic vector field or gauge choice rather than accidentally absorbing resonant `H4`.

The `11` script is the current best single-frequency quartic check.  It rebuilds `H4` directly from the operator convolution instead of using the older hand-written ordered kernel; the two differ in the truncated mode algebra, so the older ordered-kernel route should not be trusted for the Stokes test.  With direct `H4`, the non-resonant quartic solve succeeds for the convention equivalent to `K4 = H4 - 1/2 {H3,W3}`: the full coordinate contribution `1/2 {{zeta,W3},W3} + {zeta,W4}` exactly cancels the H3-only finite-depth residual and reproduces the supplied Stokes/VWA `C[3,3]`.  The finite-depth gauge check also preserves the deep-water limit `C[3,3]=3/8`.

The `12` script freezes that accepted single-frequency convention in physical units so the quartic correction can be compared directly to the finite-depth Stokes/VWA `eta33` coefficient used by the MATLAB/C++ checks.

The `13` script is the next structural diagnostic.  It asks whether the same finite-depth `W3/W4` convention that succeeds for single-frequency `C[3,3]` also satisfies the cubic and quartic homological equations as Hamiltonian identities on the truncated mode set, and whether a small symmetric broadband toy state shows any residual mismatch at the observable level.  This is the script to run when single-frequency agreement looks good but broadband wave-group results still disagree with MF12.

The `14` script is the heavier iterative search: it keeps the non-resonant quartic solve exact, but adds checkpoint logging and searches candidate sign conventions one by one so the bottleneck can be localized if the full quartic proof stalls.

The `15` script is the fast companion screen.  It keeps the single-frequency Stokes/VWA check exact, but replaces the full quartic symbolic residual by sampled Hamiltonian-residual evaluations on a fixed normal-variable state.  Its job is to tell us quickly which convention is structurally plausible before paying for a full exact residual proof.

The `13`, `14`, and `15` scripts should now be read as diagnostics and candidate-screening tools rather than the formal quartic mainline.  They establish that the accepted direct-`H4` branch is structurally plausible and isolate the successful sign convention before the final workflow is reorganized.

The new formal quartic mainline is `16` through `19`.

The `16` script freezes the accepted finite-depth `W3` basis and locks the two sign facts used downstream: the Hamiltonian convention `H3 - {H2,W3} = 0`, and the observable convention `zeta + {zeta,W3} + 1/2 {{zeta,W3},W3}` that reproduces the current H3-only single-frequency `C33`.

The `17` script makes the quartic source explicit before any homological inversion.  It prints the direct `H4` contribution, the two separate `W3`-induced quartic pieces `{H3,W3}` and `1/2 {{H2,W3},W3}`, combines them into the raw source `S4Raw`, and then verifies that the currently accepted shorthand `H4 - 1/2 {H3,W3}` is recovered only after using the frozen cubic `W3` relation.

The `18` script performs the actual quartic solve in source space.  It projects `S4H4`, `S4W3Induced`, and `S4Raw` onto their non-resonant parts, solves separately for `W4H4` and `W4W3`, defines `W4Total = W4H4 + W4W3`, and checks the three equations `{H2,W4} = -source` first numerically on sampled states and then coefficient-by-coefficient exactly.

The `19` script is the Hamiltonian-first validation pass.  It checks the frozen cubic `W3` identity, the split quartic residuals, the full coordinate map `zeta + {zeta,W3} + 1/2 {{zeta,W3},W3} + {zeta,W4Total}`, the finite-depth single-frequency `C33`, the deep-water limit `3/8`, and a small symmetric broadband toy state.  Its role is to make explicit that the new source decomposition does not replace the accepted `11/12/15` result; it clarifies where that shorthand comes from and separates Hamiltonian consistency from any remaining observable-level broadband issues.

## File Conventions

- Prefer ASCII filenames.
- Use one subfolder per topic if this directory grows.
- Keep final results reproducible from scripts whenever practical.
