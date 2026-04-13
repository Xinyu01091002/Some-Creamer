# MATLAB Prototype Workspace

This folder contains the current MATLAB prototype workspace for Creamer/MF12 comparisons. The canonical layout is now split by role so the directional and unidirectional workflows are easier to navigate.

## Structure

- `core/`
  - shared Creamer implementation and shared harmonic-separation utilities
  - canonical home of `directional_creamer_transform.m`
- `directional/`
  - directional runners and directional plotting helpers
  - use these as the main entrypoints for 2D case studies
- `unidirectional/`
  - lighter-weight 1D focused-wave workflow
  - includes `Linear_focus_envelope.m`, the 1D four-phase separator, and the live MF12 bridge
- `validation/`
  - MF12/validation loaders and live MF12 helpers
- `output/directional/`
  - directional PNG outputs
- `output/unidirectional/`
  - unidirectional PNG outputs

The older files left directly under `MATLAB/` are retained for compatibility, but the subdirectories above are the intended places to work from.

## Current Scope

- deep water only
- directional 2D case plus a lighter unidirectional 1D testbed
- second-order and exploratory third-order reconstruction
- MF12 comparison connected either through the fixed-`Akp` validation `.mat` or through live MF12 spectral calls
- a standalone non-MEX C++ backend for the expensive directional canonical-pair lambda-flow

## Important v1 Assumption

The provided directional-wave dataset gives the linear parent surface `eta(x,y)` but not the corresponding linear surface potential `phi_s(x,y)`.
The prototype therefore reconstructs the initial `phi_s` in Fourier space using a positive-frequency deep-water convention aligned with the specified propagation direction.

This remains the main closure assumption in the current implementation and should be revisited if later comparisons require a more exact phase-consistent input pair `(eta, phi_s)`.

## Current Numerical Options

- `directional_creamer_transform.m` supports three lambda-flow models:
  - `canonical_pair`
  - `legacy_zeta_only`
  - `backward_picard_316`
- `directional_creamer_transform.m` now supports optional finite depth via `cfg.depth_h`.
  - unset or `inf`: original deep-water Creamer et al. 1989 kernels
  - finite value: Wright & Creamer 1994 `theta(k)=|k| tanh(|k|h)` and finite-depth `B/D` kernels
- The active-mode truncation and lambda-step count are currently the two most important numerical knobs for third-order agreement.
- In current testing:
  - second-order `eta` and `phi` comparisons are already fairly good
  - third-order `eta33` is much more sensitive to active-mode count than to small changes in lambda-step count
  - unidirectional finite-depth `eta20/eta22` at `k_p h = 2` agrees closely with live MF12 for the baseline focused-wave case

## C++ Backend

- The C++ source lives in `../cpp/creamer_flow/`.
- MATLAB still prepares spectra, runs four-phase separation, compares against MF12, and plots.
- C++ now selects active modes and builds the canonical-pair interaction plan internally from compact binary job files, instead of requiring MATLAB to export `dest_idx`, `Bz`, `Bphi`, and `D`.
- The C++ backend supports both deep-water and finite-depth kernels through the same `cfg.depth_h` convention as the MATLAB core.
- The directional `eta33` runner can select this backend with `creamer_backend = 'cpp'`.
- Standard FFTW was not found on PATH during the current setup check; C++ therefore still receives spectral inputs from MATLAB rather than performing FFTs itself.
- A current wide-spreading stress test with `Akp=0.18`, `alpha=8`, `spread=30`, `N_lambda=6`, and `8000` active modes finished in about `204.7 s` for four phases.
- In that case, off-center `eta33` moved closer to MF12 when increasing from `6000` to `8000` modes, but centerline `eta33` slightly overshot MF12.
- On the current 16 GB workstation, treat `8000-10000` active modes as the practical range for this full-pair-kernel implementation.
- In finite-depth unidirectional tests, weak-steepness `eta33` ratios were nearly constant with `Akp` but changed strongly with `k_p h`, suggesting a structural depth dependence in the current H3-removal-only reconstruction.

## Finite-Depth Eta33 H3+H4 Check

- `unidirectional/creamer_eta33_h3h4_single_frequency_1d.m` implements the 1D single-frequency/self-interaction third-harmonic correction:
  - `Creamer(H3)` uses the finite-depth H3-only coefficient from the lambda-flow analysis.
  - `Creamer(H3+H4)` uses the Stokes/VWA coefficient recovered by absorbing the non-resonant H4 part.
- `unidirectional/run_unidirectional_finite_depth_eta33_single_frequency_case.m` is the clean regular-wave test. It prints and plots:
  - `MF12 eta33`
  - `Creamer(H3) eta33`
  - `Creamer(H3+H4) eta33`
- `unidirectional/creamer_eta33_h3h4_triad_1d.m` is the 1D broadband triad prototype:
  - It uses the finite-depth linear relation `phi_k=-i sqrt(g/theta(k)) eta_k`, `theta(k)=|k| tanh(|k|h)`.
  - It builds direct `H4` from the `N2` convolution and uses `K4=H4-1/2{H3,W3}` with only non-resonant normal monomials absorbed.
  - Its single-frequency regression at `k h=2`, `Ak=0.05` matches the Stokes/MF12 third-harmonic bin.
- `unidirectional/run_unidirectional_finite_depth_eta2_case.m` now reports:
  - `MF12 eta33`
  - `Creamer(H3) triad eta33`
  - `Creamer(H3+H4) triad eta33`
  - flow-separated and single-frequency/self-interaction eta33 as diagnostics
- `unidirectional/creamer_eta33_h3h4_triad_1d_cpp.m` calls the standalone C++ backend `cpp/creamer_flow/build/creamer_eta33_triad_1d.exe` for the same targeted eta33 calculation.
- Performance note: the MATLAB broadband triad helper now uses a targeted eta33 evaluator instead of expanding the full quartic polynomial. On the baseline `Nx=1024`, `k_p h=2` focused case, the MATLAB targeted path takes about 45 s for `max_triad_active_modes=8` and about 173 s for `12`. The C++ targeted backend takes about 0.9 s for `8`, 3.1 s for `32`, and 45 s for `60`; `100` was still too large for a 5 minute smoke-test window on this workstation. The finite-depth eta2 runner defaults to `triad_backend='cpp'` and `max_triad_active_modes=60`.
