# `PairHBNonbonded` — Custom LAMMPS Pair Style

Coarse-grained nonbonded pair potential combining hydrophobic (Morse), electrostatic (Debye–Hückel), and directional hydrogen bond interactions, with smooth mutual exclusion between intramolecular and intermolecular hydrogen bonds. HB, Morse and electrostatic terms act on **all** nonbonded pairs, both intra- and intermolecular. Implemented by Margherita Bini.

## Overview

This pair style implements the nonbonded potential used in coarse-grained (CG) molecular dynamics simulations of disordered proteins and protein aggregation. It combines three physical interaction terms: hydrophobic attraction (Morse potential), long-range electrostatics (screened Coulomb, Debye–Hückel), and directional hydrogen bonds (a geometry-dependent Gaussian well with angular selectivity).

A key feature is the smooth mutual exclusion between intramolecular and intermolecular hydrogen bonds, enforced through continuous, purely geometry-based occupancy functions to avoid force discontinuities and to avoid relying on pre-declared HB bond topology.

## Pair style syntax

```lammps
pair_style hb/nonbonded <cutoff> [kappa] [dielectric]
pair_coeff I J D0 alpha r0 lambda epsilon_hb r_hb sigma A_barrier [cutoff]
```

### Settings

| Argument     | Description                                             |
|--------------|-----------------------------------------------------------|
| `cutoff`     | Global interaction cutoff distance (Å)                  |
| `kappa`      | Debye screening parameter κ (Å⁻¹); default 0 (no screening) |
| `dielectric` | Dielectric constant ε; default 1.0                      |

### Pair coefficients

| Parameter     | Description                                                            |
|---------------|--------------------------------------------------------------------------|
| `D0`          | Morse well depth (energy units)                                        |
| `alpha`       | Morse stiffness parameter (Å⁻¹)                                        |
| `r0`          | Morse equilibrium distance (Å)                                         |
| `lambda`      | Effective charge product q_i q_j for electrostatics                   |
| `epsilon_hb`  | Hydrogen bond well depth (energy units)                                |
| `r_hb`        | Hydrogen bond equilibrium distance (Å)                                 |
| `sigma`       | Width parameter for HB Gaussian (both distance and angular terms) (Å) |
| `A_barrier`   | Height of the short-range repulsive barrier in the HB potential        |

## Potential energy

### Hydrophobic (Morse) potential

$$U_{\mathrm{hyd}}(r) = D_0 \left[ \left(e^{-\alpha(r - r_0)} - 1\right)^2 - 1 \right]$$

This produces an attractive well of depth $D_0$ at $r = r_0$.

### Electrostatic potential

$$U_{\mathrm{monopole}}(r) = \frac{q_{qrd2e}\,\lambda_{ij}}{\varepsilon \, r} \, e^{-\kappa r}$$

where $\lambda_{ij} = q_i q_j$, $\varepsilon$ is the dielectric constant, $\kappa$ is the Debye screening parameter, and $q_{qrd2e}$ is LAMMPS' internal charge-unit conversion factor.

### Hydrogen bond potential

$$U_{\mathrm{HB}}(r_{ij}, \mathbf{t}_i, \mathbf{t}_j) =
\left[
-\varepsilon_{\mathrm{hb}}
\exp\!\left(-\frac{(r_{ij} - r_{\mathrm{hb}})^2}{\sigma^2}\right)
\exp\!\left(\frac{|\hat{\mathbf{r}}_{ij} \cdot \mathbf{t}_i| - 1}{\sigma^2}\right)
\exp\!\left(\frac{|\hat{\mathbf{r}}_{ij} \cdot \mathbf{t}_j| - 1}{\sigma^2}\right)
+ A\left(\frac{r_{\mathrm{hb}}}{2^{1/6}\,r_{ij}}\right)^{12}
\right] \cdot S_{\mathrm{inter}}(i,j)$$

where $\mathbf{t}_i$, $\mathbf{t}_j$ are unit normals to the local backbone plane at *i* and *j*, and $\hat{\mathbf{r}}_{ij}$ is the unit inter-residue vector. The local backbone normal is computed from the two backbone (bond type 1) neighbors of a bead:

$$\hat{\mathbf{t}}_i = \frac{(\mathbf{r}_i - \mathbf{r}_{i-1}) \times (\mathbf{r}_{i+1} - \mathbf{r}_i)}{|(\mathbf{r}_i - \mathbf{r}_{i-1}) \times (\mathbf{r}_{i+1} - \mathbf{r}_i)|}$$

If a bead does not have exactly two backbone-bonded neighbors, a fallback normal $(0,0,1)$ is used.

The HB potential is only evaluated for pairs that are not 1-2, 1-3, or 1-4 bonded (backbone-proximal) and are not already connected by a pre-declared HB bond (bond type > 1) — the latter contribute to occupancy (see below) but are excluded here to avoid double counting with the bonded part of the potential.

The raw repulsive term $\propto (r_{\mathrm{hb}}/r)^{12}$ diverges steeply as $r \to 0$. To prevent force spikes during fast HB formation/breaking events, which could otherwise move an atom past the `comm_modify` cutoff in a single step, the barrier switches to a linear, constant-force extrapolation below $r_{\mathrm{cap}} = 0.75\, r_{\mathrm{hb}}$, matching value and slope at $r_{\mathrm{cap}}$ so the potential remains $C^1$-continuous.

> **Note:** the same soft-capping strategy, and the same $r_{\mathrm{cap}} = 0.75\,r_{hb}$ ratio, is used for the bonded HB barrier — see `README_monomer.md`.

### HB–electrostatics switching

To avoid double counting of electrostatic contributions at hydrogen bond distances, a quintic switching function $S$ smoothly suppresses the monopole term when a hydrogen bond is active:

$$S(r_{ij}, \Omega_{ij}) = 1 - \xi_{ij}^3 \left(10 - 15\,\xi_{ij} + 6\,\xi_{ij}^2\right), \qquad \xi_{ij} = \min\!\left(1,\ \frac{\left|U_{\mathrm{HB}}^{(0)}(r_{ij}, \Omega_{ij})\right|}{\varepsilon_{\mathrm{hb}}}\right)$$

where $U_{\mathrm{HB}}^{(0)}$ is the purely attractive part of $U_{\mathrm{HB}}$ (excluding the repulsive barrier). $S \to 0$ when the HB is fully formed and $S \to 1$ otherwise. This switching is only applied when `use_hb_exclusion` is enabled (default on) and $\lambda_{ij} \neq 0$.

### Total nonbonded potential

$$U = \sum_{i>j} \left[ U_{\mathrm{hyd}}(r_{ij}) + S(r_{ij}, \Omega_{ij}) \cdot U_{\mathrm{monopole}}(r_{ij}) + U_{\mathrm{HB}}(r_{ij}, \mathbf{t}_i, \mathbf{t}_j) + U_{\mathrm{rep}}^{\mathrm{destab}}(i,j) \right]$$

summed over all nonbonded pairs within cutoff, intra- and intermolecular alike. $U_{\mathrm{rep}}^{\mathrm{destab}}$ is the intramolecular destabilization penalty described below.

## Mutual exclusion between intra- and intermolecular hydrogen bonds

Each residue has at most two hydrogen bonding sites, so a maximum of two simultaneous hydrogen bonds are allowed per residue. This is enforced through two continuous, purely geometry-based occupancy functions — no dependency on pre-declared HB bond topology.

### Occupancy

For each owned atom *i*, the neighbor list is scanned once per timestep (cached, not recomputed per pair) and split by molecule membership:

$$\mathcal{O}_{\mathrm{intra}}(i) = \min\left(2,\, \sum_{k \in \mathcal{N}(i),\ \mathrm{mol}(k)=\mathrm{mol}(i)} s_{\mathrm{dist}}(d_{ik}) \cdot s_{\mathrm{ang}}(\cos\theta_{ik})\right)$$

$$\mathcal{O}_{\mathrm{inter}}(i) = \min\left(2,\, \sum_{k \in \mathcal{N}(i),\ \mathrm{mol}(k)\neq\mathrm{mol}(i)} s_{\mathrm{dist}}(d_{ik}) \cdot s_{\mathrm{ang}}(\cos\theta_{ik})\right)$$

where $\cos\theta_{ik} = |\hat{\mathbf{n}}_i \cdot \hat{\mathbf{n}}_k|$, using the same formula for both intra and inter occupancy. For $\mathcal{O}_{\mathrm{intra}}$, neighbors *k* that are 1-2, 1-3, or 1-4 bonded to *i* are excluded (backbone-proximal beads, not genuine HB partners).

Distance switching function (quintic, 5th-order Hermite):

$$s_{\mathrm{dist}}(d) =
\begin{cases}
1 & d \leq 7\,\text{Å} \\
1 - t^3(10 - 15t + 6t^2) & 7 < d < 10\,\text{Å} \\
0 & d \geq 10\,\text{Å}
\end{cases}, \quad t = \frac{d - 7}{3}$$

Angular switching function (cubic Hermite ramp on $|\cos\theta|$):

$$s_{\mathrm{ang}}(\cos\theta) =
\begin{cases}
0 & |\cos\theta| \leq 0.65 \\
u^2(3 - 2u) & 0.65 < |\cos\theta| < 0.75 \\
1 & |\cos\theta| \geq 0.75
\end{cases}, \quad u = \frac{|\cos\theta| - 0.65}{0.10}$$

The threshold at $|\cos\theta| = 0.75$ corresponds to the `hb_parallel_threshold` parameter. Occupancies are computed once per owned atom per timestep and cached; ghost atoms keep occupancy 0.

### Mutual suppression scheme

**Inter suppresses intra.** The intramolecular occupancy of each bead is fully reduced in proportion to its intermolecular occupancy, so a formed intermolecular HB is not blocked by intramolecular contacts:

$$\tilde{\mathcal{O}}_{\mathrm{intra}}(i) = \mathcal{O}_{\mathrm{intra}}(i) \cdot \left[1 - \mathrm{ramp}_1\!\left(\mathcal{O}_{\mathrm{inter}}(i)\right)\right], \qquad
\mathrm{ramp}_1(\mathcal{O}) =
\begin{cases}
0 & \mathcal{O} \leq 0 \\
\mathcal{O}^2(3 - 2\mathcal{O}) & 0 < \mathcal{O} < 1 \\
1 & \mathcal{O} \geq 1
\end{cases}$$

**Intra gates inter, partially.** For an intermolecular pair $(i,j)$, the HB potential is scaled down when either bead's corrected intramolecular occupancy is high, via a second ramp that activates only once occupancy exceeds 0.8 and caps the suppression at 50% — never fully blocking the intermolecular HB:

$$\mathrm{ramp}_2(\mathcal{O}) =
\begin{cases}
0 & \mathcal{O} \leq 0.8 \\
u^2(3 - 2u) & 0.8 < \mathcal{O} < 1 \\
1 & \mathcal{O} \geq 1
\end{cases}, \quad u = \frac{\mathcal{O} - 0.8}{0.2}$$

$$S_{\mathrm{inter}}(i,j) =
\begin{cases}
1 - 0.5\cdot\max\!\left(\mathrm{ramp}_2(\tilde{\mathcal{O}}_{\mathrm{intra}}(i)),\, \mathrm{ramp}_2(\tilde{\mathcal{O}}_{\mathrm{intra}}(j))\right) & \text{different molecules} \\
1 & \text{same molecule (no self-suppression via this gate)}
\end{cases}$$

**Soft repulsive destabilization penalty, intramolecular pairs only.** For an intramolecular pair $(i,j)$ currently in the HB-potential evaluation, if the average intermolecular engagement of *i* and *j* exceeds 0.5, a Gaussian repulsion is added directly to that pair to actively destabilize the intramolecular contact in favor of the competing intermolecular HB:

$$\text{intra\_suppress} = \tfrac{1}{2}\left[\mathrm{ramp}_1(\mathcal{O}_{\mathrm{inter}}(i)) + \mathrm{ramp}_1(\mathcal{O}_{\mathrm{inter}}(j))\right]$$

$$U_{\mathrm{rep}}^{\mathrm{destab}}(i,j) =
\begin{cases}
\varepsilon_{\mathrm{rep}} \cdot \text{intra\_suppress} \cdot \exp\!\left(-\dfrac{(r_{ij} - r_{\mathrm{rep}})^2}{\sigma_{\mathrm{rep}}^2}\right) & \text{intra\_suppress} > 0.5 \\
0 & \text{otherwise}
\end{cases}$$

with fixed parameters $\varepsilon_{\mathrm{rep}} = 0.07$, $r_{\mathrm{rep}} = 6.5\,\text{Å}$, $\sigma_{\mathrm{rep}} = 0.8\,\text{Å}$. This penalty is applied only between same-molecule pairs, and only to the pair being evaluated (not to a wider set of bond partners).

## Requirements

- LAMMPS (any recent version supporting custom pair styles)
- `atom_style` must include molecule IDs (`pair_style hb/nonbonded` requires `molecule_flag`)
- Bond types in the topology must follow the convention:
  - Type 1: backbone (harmonic) bonds — used for local normal-vector computation
  - Type > 1: pre-declared HB-type bonds (`hb6barriernew`) — used only to exclude a pair from the HB potential (to avoid double counting with the bonded term); they are not used to compute intra-/intermolecular occupancy, which is purely geometry-based (neighbor-list scan + distance/angle switching, identical treatment for intra and inter)

## Files

| File | Description |
|------|-------------|
| `pair_hb_nonbonded.cpp` | Implementation of the pair style |
| `pair_hb_nonbonded.h`   | Header file with class declaration |
