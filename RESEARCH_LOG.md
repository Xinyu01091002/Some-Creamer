# Research Log

All times below are UTC unless noted otherwise.

## 2026-03-12

### 16:21

- Confirmed the working directory was `C:\Research\Some Creamer`.
- Confirmed the directory was not yet a git repository.
- Confirmed the current source set consisted of three markdown files in [Creamer in md](/c:/Research/Some%20Creamer/Creamer%20in%20md).

### Earlier in the same session

- Read and summarized `Creamer et al. (1989)`.
- Identified the main thesis:
  - use a canonical change of variables to remove the leading cubic nonlinearity
  - keep linear time evolution in the new variables
  - reconstruct the physical surface by a nonlinear mapping
- Noted the physically important tests performed in the paper:
  - Stokes-wave shape
  - short waves riding on long waves
  - spectral modification
  - two-dimensional examples

- Read and summarized `Wright and Creamer (1994)`.
- Identified the main extension:
  - finite depth
  - slowly varying bottoms
  - slowly varying currents
- Extracted the operational prescription:
  - transform into improved linear variables locally
  - propagate with eikonal or ray methods
  - transform back to physical variables at the target depth/current state

- Read and summarized `Taylor (2022 lecture notes)`.
- Identified the main interpretive shift:
  - Creamer is used as a nonlinear geometric mapping for random seas
  - the high-wavenumber tail of the reconstructed surface is argued to tend toward a bound-wave `k^-3` form
  - cusp-like geometry is proposed as the physical-space signature of that tail

### Working conclusions as of this entry

- Best starting point for technical understanding:
  - Zakharov Hamiltonian variables
  - absence of three-wave resonance in deep-water gravity waves
  - canonical or Lie-transform elimination of `H_3`
- Best starting point for physical interpretation:
  - long-wave-induced horizontal displacement
  - short-wave packet shift, straining, and amplitude modulation
  - distinction between free-wave content and bound-wave geometry

### Open next steps

- Build a minimal formula sheet for the 1989 paper.
- Isolate the formulas needed for `\eta_t`, surface velocity, and orbital kinematics.
- Decide how to frame the "tail is bound" question:
  - frequency tail
  - wavenumber tail
  - observed spectrum versus reconstructed instantaneous geometry
