# HB6 Barrier Bond Potential — Implementation Notes

## Potential energy

Hydrogen bonds are modeled as:

$$U_{hb} = \sum_{ij} H_0 \cdot H_1 \cdot H_2 + V_{\text{barrier}}$$

where:

$$H_0 = -\epsilon \exp\!\left(-\frac{(r_{ij} - r_{hb})^2}{\sigma_r^2}\right), \quad H_1 = \exp\!\left(\frac{|\hat{r}_{ij} \cdot \mathbf{t}_i^{\text{unit}}| - 1}{\sigma_a^2}\right), \quad H_2 = \exp\!\left(\frac{|\hat{r}_{ij} \cdot \mathbf{t}_j^{\text{unit}}| - 1}{\sigma_a^2}\right)$$

Note that $\sigma_r$ (radial width, used in $H_0$) and $\sigma_a$ (angular width, used in $H_1$, $H_2$) are **independent** parameters — they no longer share a single $\sigma$.

$$V_{\text{barrier}}(r_{ij}) =
\begin{cases}
A \left(\dfrac{r_{\text{eff}}}{r_{ij}}\right)^{12}, & r_{ij} \geq r_{\text{cap}} \\[2mm]
A \left(\dfrac{r_{\text{eff}}}{r_{\text{cap}}}\right)^{12} + \left.\dfrac{\partial V_{\text{barrier}}}{\partial r_{ij}}\right|_{r_{\text{cap}}} \left(r_{ij} - r_{\text{cap}}\right), & r_{ij} < r_{\text{cap}}
\end{cases}$$

$$r_{\text{eff}} = \frac{r_{hb}}{2^{1/6}}, \qquad r_{\text{cap}} = 0.75\, r_{hb}$$

For $r_{ij} < r_{\text{cap}}$, the raw $1/r^{12}$ barrier is replaced by its linear extrapolation from $r_{\text{cap}}$, matching both value and slope at that point (keeping $V_{\text{barrier}}$ continuous and $C^1$). This is necessary because the raw form diverges as $r_{ij} \to 0$ (force $\sim 1/r^{13}$); on an **explicitly bonded** pair, an unbounded force spike can throw the bonded atom out of comm range in a single step ("Bond atoms missing"). Capping keeps the force bounded while leaving the potential unchanged in the physically relevant range $r_{ij} \geq r_{\text{cap}}$.

## Geometry

The local normal vector at site $i$ is built from the two backbone tangents:

$$\mathbf{v}_{i1} = \mathbf{r}_i - \mathbf{r}_{i-1}, \quad \mathbf{v}_{i2} = \mathbf{r}_{i+1} - \mathbf{r}_i, \quad \mathbf{t}_i = \mathbf{v}_{i1} \times \mathbf{v}_{i2}, \quad \mathbf{t}_i^{\text{unit}} = \frac{\mathbf{t}_i}{|\mathbf{t}_i|}$$

and analogously for site $j$. The bond unit vector is $\hat{r}_{ij} = (\mathbf{r}_i - \mathbf{r}_j)/r_{ij}$.

## Forces

Under the **frozen backbone approximation** ($\mathbf{t}_i^{\text{unit}}$, $\mathbf{t}_j^{\text{unit}}$ treated as constants), the force on atom $i$ is $\mathbf{F}_i = -\partial U / \partial \mathbf{r}_i$:

$$\mathbf{F}_i = -\left[\frac{\partial H_0}{\partial r_{ij}}\hat{r}_{ij}\right] H_1 H_2 - H_0\left[\frac{\partial H_1}{\partial \mathbf{r}_i}\right] H_2 - H_0 H_1\left[\frac{\partial H_2}{\partial \mathbf{r}_i}\right] - \frac{\partial V_{\text{barrier}}}{\partial r_{ij}}\hat{r}_{ij}$$

with:

$$\frac{\partial H_0}{\partial r_{ij}} = H_0 \cdot \frac{-2(r_{ij}-r_{hb})}{\sigma_r^2}$$

$$\frac{\partial H_1}{\partial \mathbf{r}_i} = \frac{H_1 \cdot \operatorname{sgn}(\hat{r}_{ij}\cdot\mathbf{t}_i^{\text{unit}})}{\sigma_a^2\, r_{ij}} \left(\mathbf{t}_i^{\text{unit}} - (\hat{r}_{ij}\cdot\mathbf{t}_i^{\text{unit}})\,\hat{r}_{ij}\right)$$

$$\frac{\partial H_2}{\partial \mathbf{r}_i} = \frac{H_2 \cdot \operatorname{sgn}(\hat{r}_{ij}\cdot\mathbf{t}_j^{\text{unit}})}{\sigma_a^2\, r_{ij}} \left(\mathbf{t}_j^{\text{unit}} - (\hat{r}_{ij}\cdot\mathbf{t}_j^{\text{unit}})\,\hat{r}_{ij}\right)$$

$$\frac{\partial V_{\text{barrier}}}{\partial r_{ij}} =
\begin{cases}
\displaystyle -\frac{12\, V_{\text{barrier}}}{r_{ij}}, & r_{ij} \geq r_{\text{cap}} \\[3mm]
\displaystyle \left.\frac{\partial V_{\text{barrier}}}{\partial r_{ij}}\right|_{r_{\text{cap}}} = -\frac{12\, V_{\text{barrier}}(r_{\text{cap}})}{r_{\text{cap}}}, & r_{ij} < r_{\text{cap}}
\end{cases}$$

i.e. below $r_{\text{cap}}$ the slope is frozen at the value the raw ($1/r^{12}$) form has exactly at $r_{\text{cap}}$, rather than continuing to diverge.

By Newton's third law, $\mathbf{F}_j = -\mathbf{F}_i$ — no separate derivative is evaluated for $j$.

> **Note:** Because backbone neighbor forces ($i\pm1$, $j\pm1$) are neglected, energy is not exactly conserved in NVE. Use NVT or Langevin dynamics.

## Parameters

| Parameter | Role |
|-----------|------|
| `epsilon` | Well depth |
| `rhb` | Equilibrium HB distance |
| `sigma_r` | Radial width (controls tolerance of $H_0$ on distance) |
| `sigma_a` | Angular width (controls tolerance of $H_1$, $H_2$ on orientation) |
| `A` | Repulsive barrier prefactor |

`bond_coeff` now takes **6** arguments (previously 5): `type epsilon rhb sigma_r sigma_a A`.

The barrier cap ratio ($r_{\text{cap}} = 0.75\, r_{hb}$) is currently hard-coded in `compute()` and `single()` as `r_cap_ratio`, not exposed as a per-type coefficient.
