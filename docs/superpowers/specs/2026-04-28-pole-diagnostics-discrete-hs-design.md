# Pole Diagnostics And Discrete HS Design

## Context

The code simulates the finite-temperature grand-canonical two-flavor Bose-Hubbard model on a triangular lattice. The current executable samples the `b` flavor explicitly and reconstructs the `c` flavor by complex conjugation. Local update is the active Monte Carlo path; global update files remain inactive.

The next development has two connected goals:

1. Add pole diagnostics to the existing continuous Hubbard-Stratonovich implementation on `main`.
2. Develop a discrete Hubbard-Stratonovich branch that reuses the same pole diagnostics so continuous and discrete auxiliary-field schemes can be compared on the same output definitions.

## Pole Diagnostics

For a fixed auxiliary-field configuration and equal-time Green matrix

```text
G = G_phi(tau),
```

the diagnostics treat `G` as a generally complex, non-Hermitian matrix. The primary pole variables are based on the complex eigenvalues `mu_a` of `G`:

```text
z_a = 1 / mu_a
d_pole = min_a |z_a|
X_pole = -log10(d_pole)
rho_G = max_a |mu_a|
```

`z_a` places the bosonic pole at the origin of the complex plane. `d_pole`, `X_pole`, and `rho_G` are equivalent scalar indicators of proximity to that pole, up to numerical roundoff. The largest singular value

```text
smax_G = max_a s_a(G)
```

is an auxiliary non-normal amplification diagnostic. It is not the pole distance and is expected to satisfy `smax_G >= rho_G` for non-normal matrices.

The determinant contribution to the sampled two-flavor weight is recorded as

```text
log_weight = log_P_HS + 2 * sum_a log(s_a(G)).
```

The singular-value form avoids overflow and underflow. Additive constants in `log_P_HS` are omitted because the diagnostic compares configurations under the same HS convention.

## Sampling Policy

Pole diagnostics are emitted once per Monte Carlo bin. They must not be folded into `Obs_equal_calc`, which samples every equal-time observation point and would make `ZGEEV` and `ZGESVD` calls far too frequent.

Each bin writes one sample using the current `Prop%Gr` at the bin-level output point. This produces a Monte Carlo time series suitable for comparing continuous and discrete auxiliary-field pole behavior.

## Output Files

The executable continues to use fixed filenames in the current run directory. Existing observable output files and their normalization are unchanged.

Add these scalar or row-wise diagnostic files:

```text
pole_z
pole_distance
pole_x
green_spectral_radius
green_smax
log_weight
```

`pole_z` writes one row per bin:

```text
Re(z_1) Im(z_1) Re(z_2) Im(z_2) ... Re(z_Ndim) Im(z_Ndim)
```

The scalar files write one real value per bin.

## Continuous HS Implementation

On `main`, continuous HS remains the active algorithm. The existing `Conf%phi_list(ns, ii, nt)` stores real continuous fields.

The HS log weight is

```text
log_P_HS = -0.5 * sum_{ns,ii,nt} phi_list(ns,ii,nt)^2.
```

This diagnostic must be added without changing the current local update proposal, acceptance ratio, Green update, or existing observable normalization.

## Discrete HS Branch

The discrete branch changes the HS field representation and local proposal while preserving the same local rank-1 Green update structure.

The discrete transformation uses four field labels:

```text
l = -2, -1, +1, +2
gamma(+/-1) = 1 + sqrt(6) / 3
gamma(+/-2) = 1 - sqrt(6) / 3
eta(+/-1) = +/- sqrt(2 * (3 - sqrt(6)))
eta(+/-2) = +/- sqrt(2 * (3 + sqrt(6)))
```

The field stored in the configuration should be an unambiguous discrete label, not just the mapped `eta` value. A helper maps the label to `eta(l)` for the exponential operator and to `log(gamma(l))` for `log_P_HS`.

For a local flip or proposal from `l_old` to `l_new`, the Metropolis ratio uses the same determinant ratio and Green rank-1 update as the continuous code, with these replacements:

```text
phi_old -> eta(l_old)
phi_new -> eta(l_new)
ratio_HS -> gamma(l_new) / gamma(l_old)
```

The discrete-branch HS log weight is

```text
log_P_HS = sum_{ns,ii,nt} log(gamma(l_{ns,ii,nt})).
```

The branch should keep fixed filename I/O and continue writing generated outputs inside run directories.

## Benchmark And Verification

Existing physics benchmark comparisons remain based on:

```text
total_NE_DQMC = mean(num_up) + mean(num_do)
total_kinetic_DQMC = mean(kinetic) * Lq
```

for live Monte Carlo runs, with the existing fixture semantics unchanged.

Pole diagnostics add validation checks but do not alter the pass/fail physics observables. Fast structural checks should confirm:

```text
pole_z exists and has one row per bin
pole_distance, pole_x, green_spectral_radius, green_smax, and log_weight exist
all scalar diagnostic values are finite where expected
pole_x is consistent with -log10(pole_distance)
green_spectral_radius is consistent with 1 / pole_distance
```

Substantive algorithm changes in the discrete HS branch must run the live DQMC benchmark suite. The optional ED recomputation target is not part of this work.

## Worktree Strategy

Implementation should proceed in two stages:

1. Add shared pole diagnostics to `main` for the existing continuous HS implementation.
2. Create an isolated branch worktree for the discrete HS implementation and reuse the same diagnostic module and output contract.

This keeps the continuous baseline available for direct comparison while isolating the algorithmic discrete-field changes.
