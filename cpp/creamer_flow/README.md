# Standalone Creamer Flow Kernel

This folder contains a standalone C++17 executable for the expensive
Creamer spectral lambda-flow. It is not a MEX file.

MATLAB remains responsible for:
- loading `.mat` data,
- building/plotting MF12 comparisons,
- phase shifts and four-phase separation,
- FFT/IFFT reconstruction.

The C++ executable reads a binary job folder, selects active modes, builds
the canonical-pair interaction plan in-process, runs fixed-step RK4 for the
lambda-flow, and writes final spectral coefficients.

Build:

```powershell
cmake -S cpp/creamer_flow -B cpp/creamer_flow/build
cmake --build cpp/creamer_flow/build --config Release
```

The current executable name is `creamer_flow_plan.exe`; this avoids
colliding with older locked `creamer_flow.exe` binaries during interactive
MATLAB sessions on Windows.

MATLAB bridge functions write:
- `meta.bin`: int64 `[ny nx n_total n_lambda_steps preserve_mean]`
- `params.bin`: double `[dx dy g energy_fraction min_active_modes max_active_modes]`
- `zeta0.bin`, `phi0.bin`: interleaved double complex spectra in MATLAB linear order

The C++ executable now builds these formerly exported plan arrays itself:
- active zero-based MATLAB-linear indices
- zero-based destination indices
- `Bz`, `Bphi`, and `D` kernel arrays

The executable writes:
- `zeta_final.bin`
- `phi_final.bin`
- `summary.txt`
