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
- The active-mode truncation and lambda-step count are currently the two most important numerical knobs for third-order agreement.
- In current testing:
  - second-order `eta` and `phi` comparisons are already fairly good
  - third-order `eta33` is much more sensitive to active-mode count than to small changes in lambda-step count

## C++ Backend

- The C++ source lives in `../cpp/creamer_flow/`.
- MATLAB still prepares spectra, runs four-phase separation, compares against MF12, and plots.
- C++ now selects active modes and builds the canonical-pair interaction plan internally from compact binary job files, instead of requiring MATLAB to export `dest_idx`, `Bz`, `Bphi`, and `D`.
- The directional `eta33` runner can select this backend with `creamer_backend = 'cpp'`.
- Standard FFTW was not found on PATH during the current setup check; C++ therefore still receives spectral inputs from MATLAB rather than performing FFTs itself.
- A current wide-spreading stress test with `Akp=0.18`, `alpha=8`, `spread=30`, `N_lambda=6`, and `8000` active modes finished in about `204.7 s` for four phases.
- In that case, off-center `eta33` moved closer to MF12 when increasing from `6000` to `8000` modes, but centerline `eta33` slightly overshot MF12.
- On the current 16 GB workstation, treat `8000-10000` active modes as the practical range for this full-pair-kernel implementation.
