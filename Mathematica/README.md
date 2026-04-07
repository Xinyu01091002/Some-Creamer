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

## File Conventions

- Prefer ASCII filenames.
- Use one subfolder per topic if this directory grows.
- Keep final results reproducible from scripts whenever practical.
