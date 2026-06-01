# Physics And Code Map

## Model

The code simulates the two-flavor Bose-Hubbard model on a triangular lattice:

```text
H = t sum_<ij> (b_i^dagger b_j + c_i^dagger c_j + h.c.)
  + U1 sum_i (n_{b,i} + n_{c,i})^2
  + U2 sum_i (n_{b,i} - n_{c,i})^2
```

The current convention is:

- `t = 1 > 0`
- `U1 < 0`
- `U2 > 0`
- flavors `b` and `c` are also called `up` and `down`

## Hubbard-Stratonovich Fields

The Trotter decomposition uses two continuous auxiliary fields:

- `phi_1` couples to `n_b + n_c` with real coefficient `sqrt(-2 * Dtau * U1)`.
- `phi_2` couples to `n_b - n_c` with imaginary coefficient `i * sqrt(2 * Dtau * U2)`.

Equivalently, the continuous-HS coupling constants are:

```text
alpha_cont(U1) = sqrt(-2 * U1 * Dtau)
alpha_cont(U2) = i * sqrt(2 * U2 * Dtau)
```

After decoupling, the two flavor Hamiltonians are complex conjugates. The code samples the `b` flavor explicitly. The `c` flavor Green matrix is reconstructed with complex conjugation.

## Green Function Convention

For a fixed auxiliary-field configuration, the equal-time Green matrix is:

```text
G_ij(tau) = < b_i b_j^dagger >
```

In code, `Prop%Gr` stores the `b` flavor matrix. `ObserEqual_mod` constructs:

```text
Grup  = Prop%Gr
Grupc = transpose(Grup) - ZKRON
Grdo  = dconjg(Prop%Gr)
Grdoc = transpose(Grdo) - ZKRON
```

## Continuous Pole Diagnostics

The continuous-HS run writes configuration-level pole diagnostics alongside the scalar observables:

| DQMC output file | Meaning |
| --- | --- |
| `pole_z` | Re/Im pairs for `z_a = 1 / mu_a(G)` |
| `pole_distance` | `min_a |z_a|` |
| `pole_x` | `-log10(pole_distance)` |
| `green_spectral_radius` | `max_a |mu_a(G)|` |
| `green_smax` | largest singular value of `G` |
| `log_weight` | `log_P_HS + 2 * sum log(s_a(G))` |

Here `mu_a(G)` are eigenvalues of the `b`-flavor Green matrix and `s_a(G)` are its singular values. For the continuous fields,

```text
log_P_HS = -0.5 * sum phi^2
```

with normalization constants omitted from `log_weight`.

These diagnostics are sampled once per bin per MPI rank. For `MPI_NP > 1`, each rank contributes one independent configuration sample per bin, and rank 0 gathers and writes `ISIZE` rows per bin in rank order. The rows are not rank averaged. Use `MPI_NP=1` when following a single Markov-chain time series or pole-spike trace.

Pole diagnostics are configuration diagnostics. They do not change the density, number, kinetic, or occupancy observable normalizations below.

## Physical Symbols To Code Variables

| Physical symbol | Code variable or location |
| --- | --- |
| `Lx`, `Ly` | `Nlx`, `Nly` in `src/calc_basic.f90`, read from `paramC_sets.txt` |
| `Lq = Lx * Ly` | `Lq` in `CalcBasic` |
| `beta` | `Beta` in `CalcBasic` |
| `Delta tau` | `Dtau = Beta / Ltrot` in `Params_set` |
| `Ltrot` | `Ltrot` in `CalcBasic`, read from `paramC_sets.txt` |
| `t` | `RT` in `CalcBasic`, currently set to `1.d0` in `Params_set` |
| `U1`, `U2` | `RU1`, `RU2` in `CalcBasic`, read from `paramC_sets.txt` |
| `mu` | `mu` in `CalcBasic`, read from `paramC_sets.txt` |
| auxiliary field flavor index | `ns = 1` for `U1`, `ns = 2` for `U2` |
| auxiliary fields | `Conf%phi_list(ns, ii, nt)` in `src/fields.f90` |
| triangular lattice nearest-neighbor bonds | `Latt%L_bonds(ii, nb)` in `src/lattice.f90` |
| space-time bonds | `Latt%LT_bonds(iit, nb)` in `src/lattice.f90` |
| `b` flavor Green matrix | `Prop%Gr` |
| `c` flavor Green matrix | `dconjg(Prop%Gr)` |
| local update shift | `shiftLoc` in `CalcBasic`, read from `paramC_sets.txt` |
| warm-up shift | `shiftWarm(1:2)` in `CalcBasic`, read from `paramC_sets.txt` |

## Observable Map

| DQMC output file | Fortran field | Meaning |
| --- | --- | --- |
| `density_up` | `Obs%density_up` | per-site density for the `b` flavor |
| `density_do` | `Obs%density_do` | per-site density for the `c` flavor |
| `num_up` | `Obs%num_up` | total `b`-flavor particle number |
| `num_do` | `Obs%num_do` | total `c`-flavor particle number |
| `kinetic` | `Obs%kinetic` | kinetic observable accumulated with division by `Lq` in `Obs_equal_calc` |
| `doubleOcc` | `Obs%doubleOcc` | per-site cross-flavor onsite density product, `<sum_i n_b,i n_c,i> / Lq` |
| `squareOcc` | `Obs%squareOcc` | per-site half normal-ordered same-flavor onsite pair, `0.5 * <sum_i,s n_s,i (n_s,i - 1)> / Lq` |
| `numsquare_up` | `Obs%numsquare_up` | total `b`-flavor number-square estimator |
| `numsquare_do` | `Obs%numsquare_do` | total `c`-flavor number-square estimator |
| `den_upup_sub11` | `Obs%den_corr_up` after Fourier transform | `b-b` density correlation for the single orbital case |
| `den_dodo_sub11` | `Obs%den_corr_do` after Fourier transform | `c-c` density correlation for the single orbital case |
| `den_updo` | `Obs%den_corr_updo` after Fourier transform | cross-flavor density correlation |

## Benchmark Normalization

The ED reference script loops over two-species fixed-`NE` blocks and reports:

- total particle number `NE = NE_b + NE_c`
- total kinetic expectation value from both flavor hopping layers
- `doubleOcc`, `squareOcc`, `numsquare_up`, and `numsquare_do`

The DQMC code writes:

- `num_up` for the `b` flavor
- `num_do` for the `c` flavor
- `kinetic` divided by `Lq`
- `doubleOcc` and `squareOcc` with the per-site normalizations in the observable map
- `numsquare_up` and `numsquare_do` as total flavor number-square estimators

Therefore the benchmark comparison should use:

```text
total_NE_DQMC = last(num_up) + last(num_do)
total_kinetic_DQMC = last(kinetic) * Lq
doubleOcc_DQMC = last(doubleOcc)
squareOcc_DQMC = last(squareOcc)
numsquare_up_DQMC = last(numsquare_up)
numsquare_do_DQMC = last(numsquare_do)
```

For real Monte Carlo runs, use sample means rather than the last bin:

```text
total_NE_DQMC = mean(num_up) + mean(num_do)
total_kinetic_DQMC = mean(kinetic) * Lq
doubleOcc_DQMC = mean(doubleOcc)
squareOcc_DQMC = mean(squareOcc)
numsquare_up_DQMC = mean(numsquare_up)
numsquare_do_DQMC = mean(numsquare_do)
```

The live DQMC benchmark estimates the uncertainty of these means by blocking the time series. The reported `stderr` is the standard error of the mean computed from block means, not the standard deviation of raw per-bin samples.

Changing this conversion is a physics-level change and must be reflected in `benchmarks/references/*.json`, `benchmarks/dqmc_references/*.json`, and `benchmarks/README.md`.
