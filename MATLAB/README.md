# MATLAB Prototype Workspace

This folder contains the current MATLAB prototype workspace for a directional Creamer transform based on the general deep-water 1989 kernels, plus MF12 comparison scripts.

## Files

- `directional_creamer_transform.m`
  - core routine
  - input: one linear parent surface `eta_lin(x,y)` plus `x_vec`, `y_vec`
  - output: one nonlinear reconstructed surface `eta_nl(x,y)` and the corresponding reconstructed `phi_s(x,y)`
- `run_directional_creamer_case.m`
  - second-order `eta20/eta22` comparison driver against MF12
- `run_directional_creamer_phi_case.m`
  - second-order `phi20/phi22` comparison driver against MF12
- `run_directional_creamer_eta33_case.m`
  - third-order `eta33` comparison driver against MF12
- `creamer_four_phase_separation.m`
  - harmonic-sector separator used to extract `eta20`, `eta22`, `eta33`, `phi20`, and `phi22`

## Current Scope

- deep water only
- directional 2D case
- second-order and exploratory third-order reconstruction
- MF12 comparison connected through the fixed-`Akp` validation `.mat`
- no performance optimization yet

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
