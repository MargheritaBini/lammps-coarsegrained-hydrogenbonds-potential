# `PairHBNonbonded` — Custom LAMMPS Pair Style

**Author:** Margherita Bini  
**Description:** Coarse-grained nonbonded pair potential combining hydrophobic (Morse), electrostatic (Debye–Hückel), and directional hydrogen bond interactions, with smooth mutual exclusion between intramolecular and intermolecular hydrogen bonds.

---

## Overview

This pair style implements the nonbonded potential used in coarse-grained (CG) molecular dynamics simulations of disordered proteins and protein aggregation. It combines three physical interaction terms:

1. **Hydrophobic attraction** — Morse potential
2. **Long-range electrostatics** — screened Coulomb (Debye–Hückel)
3. **Directional hydrogen bonds** — geometry-dependent Gaussian well with angular selectivity

A key feature is the **smooth mutual exclusion** between intramolecular and intermolecular hydrogen bonds, enforced through continuous occupancy functions to avoid force discontinuities.

---

## Pair Style Syntax

```lammps
pair_style hb/nonbonded <cutoff> [kappa] [dielectric]
pair_coeff I J D0 alpha r0 lambda epsilon_hb r_hb sigma A_barrier [cutoff]
```

### Settings

| Argument     | Description                                             |
|--------------|---------------------------------------------------------|
| `cutoff`     | Global interaction cutoff distance (Å)                  |
| `kappa`      | Debye screening parameter κ (Å⁻¹); default 0 (no screening) |
| `dielectric` | Dielectric constant ε; default 1.0                      |

### Pair Coefficients

| Parameter     | Description                                                            |
|---------------|------------------------------------------------------------------------|
| `D0`          | Morse well depth (energy units)                                        |
| `alpha`       | Morse stiffness parameter (Å⁻¹)                                        |
| `r0`          | Morse equilibrium distance (Å)                                         |
| `lambda`      | Effective charge product q_i q_j for electrostatics                   |
| `epsilon_hb`  | Hydrogen bond well depth (energy units)                                |
| `r_hb`        | Hydrogen bond equilibrium distance (Å)                                 |
| `sigma`       | Width parameter for HB Gaussian (both distance and angular terms) (Å) |
| `A_barrier`   | Height of the short-range repulsive barrier in the HB potential        |

---

## Potential Energy Terms

### 1. Hydrophobic (Morse) potential

The hydrophobic interaction between residues *i* and *j* is modeled as:

$$U_{\mathrm{hyd}}(r) = D_0 \left[ \left(e^{-\alpha(r - r_0)} - 1\right)^2 - 1 \right]$$

This produces an attractive well of depth $D_0$ at $r = r_0$.

---

### 2. Electrostatic potential

In the long-range regime ($r \gg r_{\mathrm{hb}}$), electrostatics is described by a screened Coulomb potential:

$$U_{\mathrm{monopole}}(r) = \frac{\lambda_{ij}}{\varepsilon \, r} \, e^{-\kappa r}$$

where $\lambda_{ij} = q_i q_j$ is the product of effective charges, $\varepsilon$ is the dielectric constant, and $\kappa$ is the Debye screening length.

---

### 3. Hydrogen bond potential

In the hydrogen bond regime ($r \approx r_{\mathrm{hb}}$), the dipole–dipole interaction is represented by a directional potential:

$$U_{\mathrm{HB}}(r_{ij}, \mathbf{t}_i, \mathbf{t}_j) =
\left[
-\varepsilon_{\mathrm{hb}}
\exp\!\left(-\frac{(r_{ij} - r_{\mathrm{hb}})^2}{\sigma^2}\right)
\exp\!\left(\frac{|\hat{\mathbf{r}}_{ij} \cdot \mathbf{t}_i| - 1}{\sigma^2}\right)
\exp\!\left(\frac{|\hat{\mathbf{r}}_{ij} \cdot \mathbf{t}_j| - 1}{\sigma^2}\right)
+ A\left(\frac{r_{\mathrm{hb}}}{2^{1/6}\,r_{ij}}\right)^{12}
\right]$$

where:
- $\mathbf{t}_i$, $\mathbf{t}_j$ are unit normal vectors to the local backbone plane at residues *i* and *j*
- $\hat{\mathbf{r}}_{ij}$ is the unit vector along the inter-residue direction
- $A$ is the height of the short-range repulsive barrier

The **local backbone normal** $\hat{\mathbf{t}}_i$ is computed as the normalized cross product of the two vectors connecting bead *i* to its two backbone-bonded neighbors:

$$\hat{\mathbf{t}}_i = \frac{(\mathbf{r}_i - \mathbf{r}_{i-1}) \times (\mathbf{r}_{i+1} - \mathbf{r}_i)}{|(\mathbf{r}_i - \mathbf{r}_{i-1}) \times (\mathbf{r}_{i+1} - \mathbf{r}_i)|}$$

---

### 4. HB–electrostatics switching

To avoid double counting of electrostatic contributions at hydrogen bond distances, a quintic switching function $S$ smoothly suppresses the monopole term when a hydrogen bond is active:

$$S(r_{ij}, \Omega_{ij}) = 1 - \xi_{ij}^3 \left(10 - 15\,\xi_{ij} + 6\,\xi_{ij}^2\right)$$

$$\xi_{ij} = \frac{\left|U_{\mathrm{HB}}^{(0)}(r_{ij}, \Omega_{ij})\right|}{\varepsilon_{\mathrm{hb}}}$$

where $U_{\mathrm{HB}}^{(0)}$ is the purely attractive part of $U_{\mathrm{HB}}$ (excluding the repulsive barrier). This ensures $S \to 0$ when the HB is fully formed and $S \to 1$ otherwise.

---

### 5. Total nonbonded potential

$$U = U_{\mathrm{hyd}} + \sum_{i>j} \left[ S(r_{ij}, \Omega_{ij}) \cdot U_{\mathrm{monopole}}(r_{ij}) + U_{\mathrm{HB}}(r_{ij}, \mathbf{t}_i, \mathbf{t}_j) \cdot S_{\mathrm{inter}}(i,j) \right]$$

The factor $S_{\mathrm{inter}}(i,j)$ implements the mutual exclusion between intra- and intermolecular hydrogen bonds, described in detail in the next section.

---

## Mutual Exclusion between Intra- and Intermolecular Hydrogen Bonds

Each residue has at most two hydrogen bonding sites (one donor, one acceptor), so a maximum of two simultaneous hydrogen bonds are allowed per residue. The factor $S_{\mathrm{inter}}(i,j)$ enforces this constraint continuously.

### Intramolecular occupancy

For each atom *i*, a continuous intramolecular HB occupancy is accumulated over all intramolecular HB-type bonds (bond type > 1 in the LAMMPS topology) within the same molecule:

$$\mathcal{O}_{\mathrm{intra}}(i) = \min\left(2,\, \sum_{b \in \mathcal{B}_{\mathrm{intra}}(i)} s_{\mathrm{dist}}(d_b) \cdot s_{\mathrm{ang}}(\cos\theta_b)\right)$$

where $d_b$ is the current bond distance and $\cos\theta_b = |\hat{\mathbf{n}}_i \cdot \hat{\mathbf{n}}_b|$ is the alignment between local backbone normal vectors.

**Distance switching function** (quintic, 5th-order Hermite):

$$s_{\mathrm{dist}}(d) =
\begin{cases}
1 & d \leq 7\,\text{Å} \\
1 - t^3(10 - 15t + 6t^2) & 7 < d < 10\,\text{Å} \\
0 & d \geq 10\,\text{Å}
\end{cases}, \quad t = \frac{d - 7}{3}$$

**Angular switching function** (cubic Hermite ramp on $|\cos\theta|$):

$$s_{\mathrm{ang}}(\cos\theta) =
\begin{cases}
0 & |\cos\theta| \leq 0.65 \\
u^2(3 - 2u) & 0.65 < |\cos\theta| < 0.75 \\
1 & |\cos\theta| \geq 0.75
\end{cases}, \quad u = \frac{|\cos\theta| - 0.65}{0.10}$$

The threshold at $|\cos\theta| = 0.75$ corresponds to the `hb_parallel_threshold` parameter.

---

### Intermolecular occupancy

Analogously, the intermolecular occupancy of atom *i* is computed by scanning all neighbors from the other molecule in the neighbor list:

$$\mathcal{O}_{\mathrm{inter}}(i) = \min\left(2,\, \sum_{k \in \mathcal{B}_{\mathrm{inter}}(i)} s_{\mathrm{dist}}(d_k) \cdot s_{\mathrm{ang}}(\cos\theta_k)\right)$$

using the same switching functions as above. This is geometry-based and requires no pre-declared HB bonds in the topology.

---

### Mutual suppression scheme

The two occupancies are coupled through a two-way suppression mechanism.

**Step 1 — Inter suppresses intra.**  
The intramolecular occupancy is reduced when intermolecular bonds are present, so that a formed intermolecular HB does not artificially block further intermolecular contacts:

$$\tilde{\mathcal{O}}_{\mathrm{intra}}(i) = \mathcal{O}_{\mathrm{intra}}(i) \cdot \left[1 - \mathrm{ramp}_1\!\left(\mathcal{O}_{\mathrm{inter}}(i)\right)\right]$$

where:

$$\mathrm{ramp}_1(\mathcal{O}) =
\begin{cases}
0 & \mathcal{O} \leq 0 \\
\mathcal{O}^2(3 - 2\mathcal{O}) & 0 < \mathcal{O} < 1 \\
1 & \mathcal{O} \geq 1
\end{cases}$$

**Step 2 — Soft repulsive penalty.**  
When the intermolecular occupancy is high, a Gaussian repulsive penalty is applied to the intramolecular HB bond partners of *i*, to actively disfavour the intramolecular HB:

$$U_{\mathrm{rep}}(i, b) = \varepsilon_{\mathrm{rep}} \cdot \mathrm{ramp}_1\!\left(\mathcal{O}_{\mathrm{inter}}(i)\right) \cdot \exp\!\left(-\frac{(d_b - r_{\mathrm{rep}})^2}{\sigma_{\mathrm{rep}}^2}\right)$$

with fixed parameters $\varepsilon_{\mathrm{rep}} = 0.08$, $r_{\mathrm{rep}} = 6.5\,\text{Å}$, $\sigma_{\mathrm{rep}} = 0.8\,\text{Å}$. The same penalty is applied symmetrically to the bond partners of *j*.

**Step 3 — Intra gates inter.**  
The corrected intramolecular occupancy $\tilde{\mathcal{O}}_{\mathrm{intra}}$ is used to compute the scaling factor $S_{\mathrm{inter}}$, via a second ramp that activates only once the occupancy exceeds 0.5:

$$\mathrm{ramp}_2(\mathcal{O}) =
\begin{cases}
0 & \mathcal{O} \leq 0.5 \\
u^2(3 - 2u) & 0.5 < \mathcal{O} < 1 \\
1 & \mathcal{O} \geq 1
\end{cases}, \quad u = \frac{\mathcal{O} - 0.5}{0.5}$$

The threshold at 0.5 ensures that a single weakly formed intramolecular HB does not fully suppress the intermolecular interaction, while a strongly formed intra HB suppresses it completely. The final scaling factor is:

$$S_{\mathrm{inter}}(i,j) = 1 - \max\!\left(\mathrm{ramp}_2\!\left(\tilde{\mathcal{O}}_{\mathrm{intra}}(i)\right),\, \mathrm{ramp}_2\!\left(\tilde{\mathcal{O}}_{\mathrm{intra}}(j)\right)\right)$$

This is 1 (no suppression) when neither *i* nor *j* has significant intramolecular HB occupancy, and 0 when either bead is fully engaged in intramolecular hydrogen bonds.

---

## Requirements

- LAMMPS (any recent version supporting custom pair styles)
- `atom_style` must include **molecule IDs** (`pair_style hb/nonbonded` requires `molecule_flag`)
- Bond types in the topology must follow the convention:
  - Type 1: backbone (harmonic) bonds — used for normal vector computation
  - Type > 1: HB-type bonds (`hb6barriernew`) — used for intramolecular occupancy

---

## Files

| File | Description |
|------|-------------|
| `pair_hb_nonbonded.cpp` | Implementation of the pair style |
| `pair_hb_nonbonded.h`   | Header file with class declaration |



