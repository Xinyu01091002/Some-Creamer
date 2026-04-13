# Creamer 1989 Tutorial

This document is the main teaching path for understanding and partially reproducing the deep-water construction in:

- `Creamer, Henyey, Schult, Wright (1989)`, *Improved linear representation of ocean surface waves*

The target audience is someone who already knows nonlinear water-wave theory, but does not want to bluff their way through canonical transformations or Lie transforms.

## How to use this guide

Each module has five parts:

- `Straight idea`: plain-language explanation
- `Formula spine`: the minimum equations that carry the logic
- `Mathematica task`: a local symbolic or numerical check
- `Common trap`: the misunderstanding to avoid
- `Status`: what is directly in the paper, what has been verified locally, and what still needs external confirmation

The hard rule for this project is:

- do not guess a mathematical step if it is unclear
- first check the paper
- then check local notes
- then build the smallest Mathematica test that can falsify the interpretation
- if still unclear, consult an external primary or standard source
- if a formula is simplified, rescaled, reduced to one mode, periodicized, or replaced by a toy model, write the original paper-level form first and state the reason for the simplification explicitly

## Module 0: What the paper is trying to do

### Straight idea

The paper is not trying to invent a new wave equation. It is trying to find a better set of variables for the same Hamiltonian water-wave system.

In the physical variables `(\zeta, \phi_s)`, the equations are nonlinear. The idea is to transform to new canonical variables in which the leading non-resonant nonlinearity is removed from the Hamiltonian. Then the new variables can be evolved with a simpler, almost linear Hamiltonian, and the physical nonlinearity reappears when reconstructing the physical surface.

This already explains why:

- the transformed dynamics can be linear to the retained order
- the physical surface can still display second-order bound harmonics

The harmonics are not gone. They have been moved from the evolution law into the reconstruction map.

### Formula spine

Paper locations:

- introduction and motivation: section 1
- Hamiltonian setup: section 2
- canonical and Lie transforms: section 3
- one-dimensional reconstruction: section 4

The project-level picture is:

1. start with the physical Hamiltonian
2. expand `H = H_2 + H_3 + H_4 + ...`
3. construct a near-identity canonical transform that removes `H_3`
4. truncate the transformed Hamiltonian at quadratic order
5. evolve the new variables linearly
6. reconstruct the physical surface nonlinearly

### Mathematica task

Run [Mathematica/tutorial_1989/00_notation_and_assumptions.wl](/c:/Research/Some%20Creamer/Mathematica/tutorial_1989/00_notation_and_assumptions.wl) to set the notation and operator conventions used by the rest of the scripts.
It also fixes the paper-level meaning of `\theta(k)=|k|`, the deep-water dispersion
relation `\omega_k^2 = g|k|`, and the 1D Hilbert-transform convention used later.

### Common trap

Thinking that "removing `H_3`" means "removing second-order bound harmonics from the real surface". It does not.

### Status

- `Paper directly gives`: yes
- `Mathematica verified`: not applicable at this overview level
- `External source checked`: not needed yet
- `Working hypothesis`: none

## Module 1: Hamiltonian setup and notation

### Straight idea

The physical variables are Zakharov-type surface variables:

- `\zeta(x,t)`: surface elevation
- `\phi_s(x,t)`: velocity potential evaluated at the surface

These form a canonical pair. The Hamiltonian is expanded in powers of surface slope or wave steepness. The quadratic part gives the usual linear gravity-wave dynamics. The cubic part is the leading nonlinearity that the paper wants to eliminate by a canonical transformation.

### Formula spine

Main paper equations:

- `(2.9)` full Hamiltonian in surface variables
- `(2.10)` quadratic Hamiltonian `H_2`
- `(2.11)` deep-water operator `\theta = [-(\partial_x)^2]^{1/2}`
- `(2.12)` Fourier conventions
- `(2.13)` cubic Hamiltonian `H_3`

Conceptually:

- `H_2` gives the linear dispersion and linear wave evolution
- `H_3` is the first nonlinear correction
- no three-wave resonance means `H_3` can be removed by a near-identity canonical transform
- the decomposition `H = H_2 + H_3 + H_4 + ...` comes from expanding the
  `\zeta`-dependent surface operators `\mathbf D_x, \mathbf D_z` about the flat surface
  and then collecting terms by total degree in the canonical variables
- the key degree-counting rule used later is
  `deg{F,G} = deg(F) + deg(G) - 2`,
  because a Poisson bracket differentiates each factor once and then multiplies the results

### Mathematica task

Use [Mathematica/tutorial_1989/01_hamiltonian_toy_model.wl](/c:/Research/Some%20Creamer/Mathematica/tutorial_1989/01_hamiltonian_toy_model.wl), which now stays at the paper level and explains:

- why the full Hamiltonian in `(2.9)` is nonlinear
- why it can be expanded as `H = H_2 + H_3 + H_4 + ...`
- how `H_2` and `H_3` are defined by grouping terms by total degree after expanding
  the `\zeta`-dependent operators around the flat free surface
- why `(2.10)` is the linear Hamiltonian and `(2.13)` is the first nonlinear correction

### Common trap

Confusing the operator `\theta=|k|` with the phase `kx-\omega t`. In this paper, `\theta` in section 2 is the deep-water Fourier symbol, not a travelling-wave phase.

### Status

- `Paper directly gives`: yes
- `Mathematica verified`: the script now mirrors the paper-level logic, but does not yet
  symbolically derive the D-operator expansion itself
- `External source checked`: not yet
- `Working hypothesis`: none

## Module 2: Why remove `H_3`

### Straight idea

The paper is using a normal-form idea: if the leading nonlinear interaction is non-resonant, then a canonical change of variables can absorb it into the coordinate map. That improves long-time perturbation behaviour, because naive perturbation theory often produces secular growth.

This is why the authors emphasize that canonical perturbation theory is better behaved than just adding corrections directly in time.

### Formula spine

Main paper locations:

- discussion before section 3
- denominators in `(3.4)` and `(3.5)`

Key logic:

- solve a homological equation of the form
  - "choose a generator so that `{H_2, W}` cancels `H_3`"
- this works only when the denominators do not vanish
- for deep-water gravity waves, three-wave resonances are absent
- the generator `W` must itself be cubic at this stage, because
  `deg{H_2,W} = 2 + deg(W) - 2`, so matching the cubic order of `H_3`
  forces `deg(W)=3`
- the same counting also explains why `H_3` does transform, but its change
  starts one order higher:
  `deg{H_3,W} = 3 + 3 - 2 = 4`,
  so it is omitted when the transformed Hamiltonian is kept only through cubic order

The pedagogical version is:

- transformed Hamiltonian
  - `K = H_2 + H_3 - {H_2, W} + ...`
- choose `W` so that
  - `{H_2, W} = H_3`
- then `K = H_2 + O(4)`

### Mathematica task

Run [Mathematica/tutorial_1989/02_poisson_and_lie_tools.wl](/c:/Research/Some%20Creamer/Mathematica/tutorial_1989/02_poisson_and_lie_tools.wl), which stays at the paper level and explains:

- the field-theoretic Poisson bracket `(3.6)`
- why `(\zeta,\phi_s)` are canonical variables
- how keeping only `H_2` gives the linear equations
  `\dot{\zeta}_k = |k| \phi_k`, `\dot{\phi}_k = - g \zeta_k`
- how that immediately yields the deep-water dispersion relation `\omega_k^2 = g|k|`

This module is no longer using a finite-dimensional toy Hamiltonian as its main content.

### Common trap

Thinking that absence of three-wave resonance means "there are no quadratic nonlinear effects". The statement is narrower: it only says the leading cubic Hamiltonian term is non-resonant, so it can be removed from the transformed Hamiltonian.

### Status

- `Paper directly gives`: yes, at the level of the cancellation logic
- `Mathematica verified`: the script now verifies the paper-level canonical structure and
  linear dynamics rather than a toy normal-form cancellation
- `External source checked`: still optional for deeper normal-form background
- `Working hypothesis`: none at this stage

## Module 3: Canonical-transform bridge

### Straight idea

Before reading the paper's Lie-transform section, it helps to bridge three notions:

- canonical pair
- generating function
- Poisson-bracket flow

The paper first presents a global generating functional depending on old and new variables. This is formally valid, but awkward in practice because the transformation equations are implicit. Then it switches to a Lie-transform formulation, where the change of variables is generated by a local Hamiltonian flow in an auxiliary parameter `\lambda`.

### Formula spine

Main paper equations:

- `(3.1)` generating functional `F`
- `(3.2)` and `(3.3)` implicit transformation equations
- `(3.6)` Poisson bracket
- `(3.7)` Lie-flow equation in `\lambda`
- `(3.8)` and `(3.9)` boundary conditions at `\lambda=0` and `\lambda=1`

Teaching summary:

- global form:
  - one functional `F(old, new)` defines a canonical map
- Lie form:
  - one local generator `W` defines infinitesimal canonical flow
  - integrate that flow from `\lambda=0` to `\lambda=1`
- why start with `F` at all:
  - because a generating functional is the field-theory version of the usual
    canonical generating function from mechanics
  - it defines the whole coordinate change at once through functional derivatives
  - and canonical structure is built in from the start
- why not stop there:
  - because the resulting relations are implicit and awkward for order-by-order work
  - Lie transform rewrites the same canonical idea in a differential form that is
    better suited to cancelling `H_3`
- what to keep in mind about Lie transform:
  - it is still a calculational tool, not a new physical theory
  - it is preferred here because the flow equation and Poisson-bracket counting
    make the cubic cancellation logic transparent
- why `(3.1)` does not include every cubic combination:
  - the paper explicitly says it is not the most general cubic form
  - other combinations would generate structures not present in the original Hamiltonian
  - the authors keep the most general cubic form odd in `\phi_s`, so that the transformed
    Hamiltonian remains even in `\phi_s`
- why the Lie-transform section reuses the same `B,D`:
  - because, to the first non-trivial order, the Lie-transform relations are the same as
    those from the global generating functional when `W` is taken to be the negative of
    the cubic part of `(3.1)`
  - so the kernels already solved from the global method can be reused at cubic order
  - the real difference between the two schemes shows up in how higher-order structure is
    incorporated, not in the existence of a completely different third-order kernel set

### Mathematica task

Run [Mathematica/tutorial_1989/03_global_to_lie_bridge.wl](/c:/Research/Some%20Creamer/Mathematica/tutorial_1989/03_global_to_lie_bridge.wl), which now summarizes:

- the exact paper-level generating functional `(3.1)`
- the transformation equations `(3.2)-(3.3)`
- the Lie-flow definition `(3.7)`
- the cancellation logic `(3.17)-(3.19)`

The script now stays with the paper's own symbols and formulas instead of using a
one-degree-of-freedom canonical analogue.

### Common trap

Treating the parameter `\lambda` as physical time. In the paper it is an auxiliary deformation parameter, not the wave-evolution time.

### Status

- `Paper directly gives`: yes
- `Mathematica verified`: the script now restates the exact paper-level objects and
  their logical roles, rather than using a toy canonical template
- `External source checked`: still useful for general Lie-transform language if needed
- `Working hypothesis`: none at this stage

## Module 4: From global generating function to Lie transform

### Straight idea

The paper's main practical move is:

- keep the canonical-transform idea
- replace the clumsy global map with a sequence of infinitesimal local transforms

This makes the transformation easier to iterate and easier to connect to one-dimensional geometry later.

The key cancellation appears explicitly in the transformed Hamiltonian:

- `K = H_2 + H_3 - {H_2, W} + ...`

With the right `W`, the last two terms cancel.

### Formula spine

Main paper equations:

- `(3.11)` Lie generator `W`
- `(3.12)` evolution equations for `Z` and `\Phi`
- `(3.15)` and `(3.16)` transformed-Hamiltonian logic
- `(3.17)` transformed Hamiltonian through third order
- `(3.18)` explicit cancellation
- `(3.19)` transformed Hamiltonian truncated to quadratic order

What matters most conceptually is:

- the transformed Hamiltonian is simpler
- the transformed variables are still canonical
- the reconstruction from linear variables back to physical variables is where physical bound structure reappears

### Mathematica task

Run [Mathematica/tutorial_1989/03_global_to_lie_bridge.wl](/c:/Research/Some%20Creamer/Mathematica/tutorial_1989/03_global_to_lie_bridge.wl) as a paper-level summary of:

- the generating functional `(3.1)`
- the transformation equations `(3.2)-(3.3)`
- the Lie-generator `W` and lambda-flow `(3.7)`
- the cubic-term cancellation `(3.17)-(3.19)`

The point is not to rederive `(3.4)` and `(3.5)` symbolically. The point is to keep
the exact paper-level logic visible before any later experiment or simplification.

### Common trap

Expecting a full symbolic re-derivation of the paper's kernels `B(1,2,3)` and `D(1,2,3)` from Mathematica alone. That is not the goal of the first stage.

### Status

- `Paper directly gives`: yes
- `Mathematica verified`: the script now mirrors the paper-level cancellation logic
  without using a reduced canonical toy model
- `External source checked`: recommended if the abstract Lie-series machinery still feels opaque
- `Working hypothesis`: none at this stage

## Module 5: One-dimensional specialization and Hilbert geometry

### Straight idea

Deep-water one dimension is special. The Lie-transform equations simplify so much that the transformed field can be interpreted in terms of the Hilbert transform of the surface profile. That Hilbert transform plays the role of an approximate horizontal Lagrangian displacement.

This is where the geometric picture becomes concrete:

- the new variables evolve linearly
- the physical surface is obtained by a horizontal remapping
- this remapping generates bound harmonics and crest sharpening

### Formula spine

Main paper equations:

- `(4.1)` one-dimensional simplification of `B`
- `(4.2)` and `(4.3)` Hilbert transform definition
- `(4.7)` one-dimensional Lie-flow equations
- `(4.8)` and `(4.9)` characteristic solution and shifted coordinate
- `(4.10)` inverse-direction reconstruction form

The essential interpretation is:

- `\tilde Z` is approximately the horizontal displacement
- `x` and `\chi` differ by that displacement
- a simple profile in the shifted coordinate becomes a nonlinear Eulerian profile after remapping
- the algebraic reason the 1D case is special is that `1D + k_1+k_2+k_3=0` already forces
  one `|k|` to equal the sum of the other two, so the general `B,D` kernels collapse
- this is a consequence of one-dimensional triad closure itself, not a separately imposed
  ordering of `k_1, k_2, k_3`
- `\chi` should be read as the characteristic label of a surface point, so section 4 is
  really solving the coordinate remapping hidden inside the abstract canonical transform

### Mathematica task

Run [Mathematica/tutorial_1989/04_1d_hilbert_remap.wl](/c:/Research/Some%20Creamer/Mathematica/tutorial_1989/04_1d_hilbert_remap.wl). It checks two things:

- the Hilbert-transform symbol in Fourier space
- a remapping-generated harmonic expansion for a single travelling wave

### Common trap

Thinking the Hilbert-transform picture is a generic fact about all depths or all canonical transforms. The paper's own 1994 finite-depth extension strongly suggests that this simple geometry is a deep-water special simplification.

### Status

- `Paper directly gives`: yes
- `Mathematica verified`: yes, for the harmonic-generation picture and Fourier-sign convention
- `External source checked`: not required yet
- `Working hypothesis`: the remapping interpretation is exact in the 1D deep-water Lie-transform framework used here, but should not be blindly exported to finite depth

## Module 6: Reconstructing the physical surface to `(4.13)`

### Straight idea

This is the payoff of stage 1. The paper shows how to go from the linear variables of the truncated Hamiltonian back to the physical surface. The result is the nonlinear reconstruction formula `(4.13)`.

That formula is the cleanest place to understand the Creamer transform in practice:

- linear evolution lives in the transformed variables
- the physical surface is nonlinear because the map back is nonlinear
- bound structure is generated during reconstruction, not as a separate free-wave evolution

### Formula spine

Main paper equations:

- `(4.11)` and `(4.12)` linear variables in terms of physical variables
- `(4.13)` physical variables in terms of the linear variables
- `(4.14)` Fourier-space reconstruction
- `(4.15)` derivatives
- `(4.18)` and `(4.19)` velocity reconstruction ingredients

For stage 1, the main conceptual target is:

- understand `(4.13)` and `(4.14)` as reconstruction formulas
- understand why they naturally create bound harmonics
- understand `(4.13)` as the final 1D expression of the already-solved coordinate remapping:
  linear transformed variables are pushed back to the physical surface by the `x <-> chi`
  characteristic map
- understand that "reconstruction" here does not mean solving a new dynamics:
  it means translating the same surface-point set from parameter coordinates back into
  Eulerian surface coordinates

### Mathematica task

Run [Mathematica/tutorial_1989/05_reconstruct_4_13.wl](/c:/Research/Some%20Creamer/Mathematica/tutorial_1989/05_reconstruct_4_13.wl). It uses a finite-period pedagogical adaptation of the Fourier reconstruction idea behind `(4.14a)` to compare:

- a linear input profile
- a reconstructed nonlinear profile generated from the same parent field

The finite-period normalization in the script is a teaching adaptation. It is not claimed to be the paper's exact continuum normalization.

### Common trap

Confusing "reconstructed profile contains high harmonics" with "the model introduced new independent free waves". In this stage, the extra harmonics are reconstruction-generated bound structure.

### Status

- `Paper directly gives`: yes
- `Mathematica verified`: yes, in a periodic teaching adaptation of the reconstruction formula
- `External source checked`: not required yet
- `Working hypothesis`: the script captures the reconstruction mechanism while using periodic rather than continuum normalization

## Stage-1 outcome

At the end of stage 1, the intended understanding is:

- `Creamer transform` is best viewed as a near-identity canonical change of variables on Zakharov-type free-surface variables
- its main Hamiltonian purpose is to remove the cubic non-resonant term `H_3`
- in deep-water 1D, the Lie-transform version becomes geometrically intelligible through a Hilbert-transform-based horizontal displacement
- the physical wave shape is then recovered by a nonlinear reconstruction map
- the extra harmonics seen in the physical surface are therefore mainly bound structure created by reconstruction

## What remains for later stages

- reproduce the single-frequency Stokes comparison more quantitatively
- reproduce the short-wave packet on long-wave analysis from section 5
- examine how much of the deep-water geometry survives in finite depth
- decide whether the finite-depth "best mapping" problem should be attacked by higher-order normal forms, better reconstruction maps, or both
