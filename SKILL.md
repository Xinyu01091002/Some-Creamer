# Creamer Transform Working Guide

## Purpose

Use this file as a compact project memory for future work in this directory.

The working objective is not just to read the literature, but to turn the Creamer transform into a usable mental and computational tool for:

- surface reconstruction
- kinematics
- time derivatives and surface velocities
- spectral interpretation
- the bound-versus-free status of tails in frequency and wavenumber space

## Core Idea

The Creamer transform is a near-identity canonical transformation of the free-surface variables, usually written in terms of Zakharov variables `(\zeta, \phi_s)`.

The point of the transformation is:

- remove the cubic Hamiltonian term when three-wave resonances are absent
- evolve the transformed variables with linear dynamics
- recover much of the leading nonlinear physics when mapping back to the physical surface

In 1-D, the Lie-transform version admits a useful kinematic interpretation:

- the Hilbert transform of the linear surface is approximately the horizontal Lagrangian displacement
- the nonlinear surface is obtained by a horizontal remapping plus amplitude reorganization

## What Matters Most

### 1989 paper

- Deep-water formulation.
- Canonical/Lie-transform removal of `H_3`.
- Strong practical emphasis on:
  - Stokes-wave approximation
  - short-wave packets on long waves
  - spectral correction
  - 2-D behavior

### 1994 paper

- Finite-depth extension.
- Replace deep-water operator by finite-depth `\theta(k) = k \tanh(kh)`.
- Combine local nonlinear transformation with linear eikonal transport over slowly varying bathymetry or currents.

### Taylor notes

- Recast Creamer as a nonlinear geometry map for random seas.
- Emphasize that high-`k` tails in the reconstructed surface may be predominantly bound.
- Connect the universal `k^-3` tail to cusp-like local geometry.

## Minimal Knowledge Prerequisites

Assume the reader already knows water-wave theory. Then the needed extras are:

1. Zakharov Hamiltonian formulation for free-surface gravity waves.
2. Canonical perturbation theory and Lie transforms.
3. Why absence of 3-wave resonance allows elimination of `H_3`.
4. Difference between free and bound waves in a weakly nonlinear random sea.
5. How long-wave-induced advection and straining act on short waves.
6. How spectral tails are connected to local regularity or singularity of the surface profile.

## Working Heuristics

- Think of the transform as a better coordinate system, not as a new dynamical law.
- Separate "linear evolution in transformed variables" from "nonlinear reconstruction in physical variables".
- For short-on-long interaction, expect the transform to encode:
  - packet displacement
  - wavelength modulation
  - envelope compression or stretching
- For spectral questions, distinguish carefully:
  - the input linear spectrum
  - the reconstructed physical-surface spectrum
  - free-wave transport versus bound-wave geometry

## Next Technical Targets

1. Build a one-page formula sheet from the 1989 paper.
2. Write out the cleanest route from transformed variables to `\eta_t`.
3. Derive surface horizontal and vertical velocities in the reconstructed field.
4. Compare the meaning of tails in:
   - `S(\omega)`
   - `S(k)`
   - instantaneous reconstructed profiles
5. Test whether the tail question should be asked in terms of:
   - propagation
   - phase locking
   - canonical normal form
   - observational diagnostics

## Repository Intent

This folder should evolve toward a small research workspace that supports:

- reading notes
- formula extraction
- derivations
- numerical experiments
- project-level conclusions about possible applications of the Creamer transform
