# CLAUDE.md — AI Assistant Guide for Some-Creamer

This file provides context for AI assistants working in this repository. It covers project purpose, codebase structure, development workflows, and conventions to follow.

---

## Project Overview

**Some-Creamer** is an active research workspace for understanding and numerically implementing the **Creamer canonical transformation** for water-wave theory. The core scientific goal is to:

1. Understand the Creamer transform as a near-identity canonical transformation on Zakharov variables `(ζ, φ_s)` that removes the cubic Hamiltonian term `H_3` when three-wave resonances are absent.
2. Implement the transform numerically to reconstruct nonlinear wave surfaces from linear parent fields.
3. Validate against MF12 (Maltby-Folley 2012 Hamiltonian wave model) predictions.
4. Extend from deep-water (Creamer et al. 1989) to finite-depth regimes (Wright & Creamer 1994).

This is **not** a software product — it is a research codebase. Expect research-grade code and prioritize correctness and scientific clarity over software engineering polish.

---

## Repository Structure

```
Some-Creamer/
├── CLAUDE.md                            # This file
├── README.md                            # Project overview and current status
├── SKILL.md                             # Compact working guide and heuristics
├── RESEARCH_LOG.md                      # Dated session log (2026-03-12 onward)
├── Creamer_1989_Tutorial.md             # Staged tutorial for the 1989 paper
├── Foundations_Canonical_Variables.md   # Bridge note on canonical variables
├── Lie_Transform_Notes.md               # Bridge note on Lie transforms
├── Practical_1D_Creamer_Form.md         # 1D kinematic interpretation
├── Creamer in md/                       # Markdown-converted source papers
│   ├── Creamer et al. (1989)
│   ├── Wright and Creamer (1994)
│   └── Taylor (2022 notes)
├── MATLAB/                              # Primary numerical workspace
│   ├── core/                            # Shared Creamer implementation
│   ├── directional/                     # 2D directional case runners
│   ├── unidirectional/                  # 1D focused-wave testbed
│   ├── validation/                      # MF12 validation bridges
│   └── output/                          # Plot output (excluded from git)
├── Mathematica/                         # Symbolic derivations
│   ├── tutorial_1989/                   # 6 scripts for the 1989 paper
│   └── finite_depth/                    # 19 scripts for finite-depth theory
└── cpp/creamer_flow/                    # C++17 performance backend
    ├── CMakeLists.txt
    ├── creamer_flow.cpp                 # Main lambda-flow executable
    └── creamer_eta33_triad_1d.cpp       # 1D triad analysis
```

---

## Architecture: Three-Tier Execution

The numerical workflow is split across three tiers:

| Tier | Technology | Role |
|------|-----------|------|
| Orchestration | MATLAB | Data I/O, FFT/IFFT, MF12 comparison, four-phase separation, plotting |
| Symbolic | Wolfram Mathematica (`.wl` scripts) | Perturbation series, operator algebra, Lie-transform derivations |
| Compute-intensive | C++17 | Canonical-pair lambda-flow (RK4), active-mode interaction plan, OpenMP |

**Data flow summary:**
1. Load linear parent surface `η_lin` (MATLAB)
2. FFT → Fourier space `ζ̂`, reconstruct `φ̄̂` using depth-dependent `θ(k)`
3. Select active modes by energy fraction
4. Build interaction plan (B, D kernels, triple indices)
5. Lambda-flow RK4 integration: λ: 1 → 0 (C++ or MATLAB)
6. Inverse map back to physical space (IFFT)
7. Four-phase separation: extract `η20, η22, η33`, `φ20, φ22`
8. Compare against MF12 reference

---

## MATLAB Codebase

### Core modules (`MATLAB/core/`)

| File | Purpose |
|------|---------|
| `directional_creamer_transform.m` | Main 2D Creamer transform; supports `canonical_pair`, `legacy_zeta_only`, `backward_picard_316` lambda-flow variants |
| `directional_creamer_transform_cpp.m` | Wrapper that delegates to C++ backend |
| `creamer_four_phase_separation.m` | Harmonic decomposition; separates `η20, η22, η33, φ20, φ22` |
| `creamer_four_phase_separation_cpp.m` | C++ version of four-phase separator |
| `directional_frequency_filter.m` | Frequency filtering utilities |
| `apply_global_phase_shift_2d.m` | Phase correction |

### Runner scripts

**Directional (2D) runners** (`MATLAB/directional/`):
- `run_directional_creamer_case.m` — locked case `Akp=0.12`, `alpha=8`, `chi=0°`, `spread=30°`; validates `η20/η22`
- `run_directional_creamer_eta33_case.m` — third-order harmonic analysis
- `run_directional_creamer_phi_case.m` — surface potential `φ20/φ22` reconstruction
- `plot_creamer_centerline_spectrum.m` — spectral visualization

**Unidirectional (1D) runners** (`MATLAB/unidirectional/`):
- `run_unidirectional_creamer_case.m` — Gaussian focused-wave testbed; live MF12 comparison
- `run_unidirectional_finite_depth_eta2_case.m` — finite-depth baseline at `k_p h=2`
- `run_unidirectional_finite_depth_eta33_single_frequency_case.m` — single-frequency check
- `creamer_eta33_h3h4_single_frequency_1d.m` — H3+H4 correction prototype
- `creamer_eta33_h3h4_triad_1d.m` — broadband triad analysis
- `creamer_eta33_h3h4_triad_1d_cpp.m` — C++ backend for triad

**Validation** (`MATLAB/validation/`):
- `load_validation_second_order_case.m`, `load_validation_third_order_case.m`
- `match_validation_second_order_case.m`
- `mf12_from_linear_surface.m` — MF12 computation interface

### MATLAB Conventions

- **Naming**: `snake_case` for all function and variable names.
- **Mathematical names**: Greek symbols spelled out (`zeta`, `phi_s`, `eta`, `theta_k`).
- **Private helpers**: prefixed with `local_` inside a file.
- **Configuration**: passed as `cfg.*` struct — e.g. `cfg.depth_h`, `cfg.lambda_flow_variant`.
- **Depth mode**: `cfg.depth_h = inf` → deep-water kernels; finite value → Wright & Creamer 1994 kernels with `θ(k) = |k| tanh(|k|h)`.
- **Physical constants**: `g = 9.81` m/s².
- **No build step** — MATLAB files run directly.

---

## C++ Backend (`cpp/creamer_flow/`)

### Build

```bash
cmake -S cpp/creamer_flow -B cpp/creamer_flow/build
cmake --build cpp/creamer_flow/build --config Release
```

Produces two executables:
- `creamer_flow_plan` — main lambda-flow binary
- `creamer_eta33_triad_1d` — standalone 1D triad analyzer

**Requirements:** C++17 compiler, CMake 3.16+. OpenMP optional (auto-detected).
**Optimization:** `-O3 -march=native` on GCC/Clang; `/O2` on MSVC.
**Excluded from git:** `cpp/**/build/`

### C++ I/O Protocol

The main executable uses binary file I/O (not MEX):

| File | Direction | Content |
|------|-----------|---------|
| `meta.bin` | Input | Grid/job metadata |
| `params.bin` | Input | Physical parameters |
| `zeta0.bin` | Input | Initial ζ field |
| `phi0.bin` | Input | Initial φ field |
| `zeta_final.bin` | Output | Final ζ after lambda-flow |
| `phi_final.bin` | Output | Final φ after lambda-flow |
| `summary.txt` | Output | Run summary / diagnostics |

### C++ Conventions

- **Standard**: C++17, no external numerical libraries — all kernels are implemented inline.
- **Parallelism**: OpenMP `#pragma omp parallel for` over active-mode loops.
- **Active-mode limit**: Practical range is 8000–10000 modes on a 16 GB workstation. Beyond that, switch to block/on-the-fly interaction kernels rather than storing all pair arrays.
- **Performance reference**: Triad analyzer: ~0.9 s for 8 active modes, ~45 s for 60 active modes.

---

## Mathematica Symbolic Scripts

### Conventions

- All scripts are `.wl` files runnable with `wolframscript`.
- **Tutorial series** (`tutorial_1989/`, numbered `00`–`05`): work through in order; covers notation, toy Hamiltonian, Poisson tools, Lie bridge, Hilbert remap, and surface reconstruction.
- **Finite-depth series** (`finite_depth/`, numbered `01`–`19`): staged investigation of H3+H4 absorption for finite-depth regimes.
- Prefer `.wl` scripts for successful derivations; use `.nb` notebooks for exploratory work.
- Keep symbolic outputs synced back into the relevant markdown notes.
- **Strategy**: Break hard derivations into small tractable symbolic tasks — never try to automate the full Lie-transform engine; keep theory logic explicit.

---

## Documentation Conventions

The repository maintains several layered documentation files with distinct roles:

| File | Role |
|------|------|
| `README.md` | Current project status, research questions, MATLAB status |
| `SKILL.md` | Compact AI working guide — heuristics and conceptual framework |
| `RESEARCH_LOG.md` | Dated session log — always append, never rewrite history |
| `Creamer_1989_Tutorial.md` | Pedagogical walkthrough of the 1989 paper |
| `Foundations_Canonical_Variables.md` | Bridge note: canonical variables for AI/reader onboarding |
| `Lie_Transform_Notes.md` | Bridge note: Lie transforms and λ-flow |
| `Practical_1D_Creamer_Form.md` | 1D kinematic interpretation |

**When to update documentation:**
- `RESEARCH_LOG.md`: After each session or significant result. Append with a dated section.
- `README.md`: When the current status, MATLAB structure, or open questions change materially.
- `SKILL.md`: When a new working heuristic or conceptual clarification has been validated.
- Bridge notes: When a conceptual gap is resolved that would have blocked earlier work.

---

## Development Workflow

### Typical Session Pattern

1. Read `SKILL.md` to orient to current heuristics and open questions.
2. Check `RESEARCH_LOG.md` for the most recent session state.
3. Run the relevant MATLAB runner from `directional/` or `unidirectional/`.
4. If compute-heavy, build and invoke the C++ backend; results feed back into MATLAB for plotting.
5. Use Mathematica for symbolic work on new theory targets.
6. Append findings to `RESEARCH_LOG.md` and update `README.md` if MATLAB status or open questions change.

### Git Conventions

- **Branch**: development work goes on feature branches; `master` holds stable milestones.
- **Current active branch**: `claude/add-claude-documentation-YJKfU`
- **Commit messages**: descriptive imperative style — e.g. `Add H3+H4 eta33 prototype, finite-depth Mathematica scripts`.
- **Excluded files** (`.gitignore`): `*.mat`, `*.png`, `*.asv`, `*.m~`, `tmp_*`, `external/`, `build/`
- **Large datasets**: never commit `.mat` files — document their filenames and generation scripts in `README.md` instead.

### Validation Approach

No automated CI/CD exists. Validation is manual:

| Test | How |
|------|-----|
| 2nd-order `η20, η22` | Compare Creamer output vs. MF12 from `.mat` dataset; should agree well |
| 2nd-order `φ20, φ22` | Same dataset; currently "quite good" agreement |
| 3rd-order `η33` directional | Sensitive to active-mode count; agree well for narrow-spreading cases |
| Unidirectional 1D | Focused Gaussian packets; live MF12 comparison at multiple `Akp` |
| Finite-depth | `k_p h=2` baseline; expect systematic depth-dependent gap in `η33` |

**Known systematic gap:** For finite depth, `Creamer η33 / MF12 η33 ≈ 0.10` at `k_p h=1` and `≈ 0.19` at `k_p h=2`. This is a depth-structure issue (H3+H4 absorption), not a steepness or active-mode issue.

---

## Key Physical Concepts for AI Reasoning

When reasoning about this codebase, keep these conceptual distinctions in mind:

1. **Canonical vs. physical variables**: The transform works in Zakharov variables `(ζ, φ_s)`; the physical surface `η` is recovered via a nonlinear map.

2. **H3 removal ≠ removing second-order physics**: Eliminating `H_3` from the transformed Hamiltonian does not remove second-order bound harmonics from the physical surface. It relocates their origin from the evolution equations into the reconstruction map.

3. **Deep vs. finite depth**:
   - Deep water: H3 is globally non-resonant → clean removal; simple Hilbert-transform geometric remapping picture.
   - Finite depth: H3 contains resonant, near-resonant, and non-resonant parts. Only the non-resonant part is cleanly removable; higher-order Hamiltonian terms (H4) may be essential.

4. **λ-flow**: The canonical transformation is implemented as a parameter-flow λ: 1→0, integrated with RK4. This is not time evolution — it is a deformation of variables.

5. **Active-mode selection**: Before building the interaction plan, modes are selected by an energy-fraction threshold. Too few modes → poor η33; too many → prohibitive compute cost.

6. **Four-phase separation**: After reconstruction, harmonic components are extracted by phase:
   - `η20`: zeroth-harmonic second-order (set-down/set-up)
   - `η22`: second-harmonic second-order (Stokes-like)
   - `η33`: third-harmonic third-order
   - `φ20`, `φ22`: corresponding velocity potential harmonics

---

## What NOT to Do

- **Do not add engineering abstractions** not required by the physics (factories, dependency injection, etc.).
- **Do not commit `.mat` datasets** — document them in README and generate them with runner scripts.
- **Do not rewrite RESEARCH_LOG.md** — it is a history log; always append.
- **Do not change the λ-flow integration scheme** (RK4) without explicitly documenting the theoretical justification in the log.
- **Do not assume deep-water results transfer to finite depth** — always check whether the geometric remapping picture still holds.
- **Do not skip the four-phase separation step** when comparing with MF12 — MF12 provides harmonic-decomposed outputs.

---

## Quick Reference: Running the Code

### MATLAB (directional second-order)
```matlab
% From MATLAB/directional/
run_directional_creamer_case
```

### MATLAB (unidirectional, live MF12)
```matlab
% From MATLAB/unidirectional/
run_unidirectional_creamer_case
```

### C++ backend: build and run
```bash
# Build
cmake -S cpp/creamer_flow -B cpp/creamer_flow/build
cmake --build cpp/creamer_flow/build --config Release

# Run (called from MATLAB via directional_creamer_transform_cpp.m)
# Executable: cpp/creamer_flow/build/creamer_flow_plan
```

### Mathematica symbolic derivations
```bash
# Run a script from terminal
wolframscript -file Mathematica/tutorial_1989/01_hamiltonian_toy_model.wl
```

---

## Current Status (as of 2026-04-13)

| Component | Status |
|-----------|--------|
| Deep-water η20/η22 vs MF12 | Validated — "quite good" agreement |
| Deep-water φ20/φ22 vs MF12 | Validated — "quite good" agreement |
| Directional η33 (narrow spread) | Good when enough active modes are used |
| Directional η33 (wide spread) | Still sensitive to active-mode count (~8000–10000 needed) |
| Finite-depth η20/η22 | Validated at k_p h=2 |
| Finite-depth η33 | Systematic gap (~10–19% of MF12) — H3+H4 absorption in progress |
| Mathematica H4 derivation | 19 scripts underway (finite_depth/01–19) |
| C++ backend | Stable for canonical_pair lambda-flow + triad analysis |

**Active research front:** Finite-depth H3+H4 normal-form closure — whether removing only H3 is sufficient or whether H4 quartic absorption is needed to close the η33 gap.
