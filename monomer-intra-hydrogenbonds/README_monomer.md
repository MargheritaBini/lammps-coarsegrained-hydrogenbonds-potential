# HB6 Barrier Bond Potential — Implementation Notes

## Potential energy

Hydrogen bonds are modeled as:

$$U_{hb} = \sum_{ij} H_0 \cdot H_1 \cdot H_2 + V_{\text{barrier}}$$

where:

$$H_0 = -\epsilon \exp\!\left(-\frac{(r_{ij} - r_{hb})^2}{\sigma^2}\right), \quad H_1 = \exp\!\left(\frac{|\hat{r}_{ij} \cdot \mathbf{t}_i^{\text{unit}}| - 1}{\sigma^2}\right), \quad H_2 = \exp\!\left(\frac{|\hat{r}_{ij} \cdot \mathbf{t}_j^{\text{unit}}| - 1}{\sigma^2}\right)$$

$$V_{\text{barrier}} = A \left(\frac{r_{\text{eff}}}{r_{ij}}\right)^{12}, \qquad r_{\text{eff}} = \frac{r_{hb}}{2^{1/6}}$$

## Geometry

The local normal vector at site $i$ is built from the two backbone tangents:

$$\mathbf{v}_{i1} = \mathbf{r}_i - \mathbf{r}_{i-1}, \quad \mathbf{v}_{i2} = \mathbf{r}_{i+1} - \mathbf{r}_i, \quad \mathbf{t}_i = \mathbf{v}_{i1} \times \mathbf{v}_{i2}, \quad \mathbf{t}_i^{\text{unit}} = \frac{\mathbf{t}_i}{|\mathbf{t}_i|}$$

and analogously for site $j$. The bond unit vector is $\hat{r}_{ij} = (\mathbf{r}_i - \mathbf{r}_j)/r_{ij}$.

## Forces

Under the **frozen backbone approximation** ($\mathbf{t}_i^{\text{unit}}$, $\mathbf{t}_j^{\text{unit}}$ treated as constants), the force on atom $i$ is $\mathbf{F}_i = -\partial U / \partial \mathbf{r}_i$:

$$\mathbf{F}_i = -\left[\frac{\partial H_0}{\partial r_{ij}}\hat{r}_{ij}\right] H_1 H_2 - H_0\left[\frac{\partial H_1}{\partial \mathbf{r}_i}\right] H_2 - H_0 H_1\left[\frac{\partial H_2}{\partial \mathbf{r}_i}\right] - \frac{\partial V_{\text{barrier}}}{\partial r_{ij}}\hat{r}_{ij}$$

with:

$$\frac{\partial H_0}{\partial r_{ij}} = H_0 \cdot \frac{-2(r_{ij}-r_{hb})}{\sigma^2}$$

$$\frac{\partial H_1}{\partial \mathbf{r}_i} = \frac{H_1 \cdot \operatorname{sgn}(\hat{r}_{ij}\cdot\mathbf{t}_i^{\text{unit}})}{\sigma^2\, r_{ij}} \left(\mathbf{t}_i^{\text{unit}} - (\hat{r}_{ij}\cdot\mathbf{t}_i^{\text{unit}})\,\hat{r}_{ij}\right)$$

$$\frac{\partial H_2}{\partial \mathbf{r}_i} = \frac{H_2 \cdot \operatorname{sgn}(\hat{r}_{ij}\cdot\mathbf{t}_j^{\text{unit}})}{\sigma^2\, r_{ij}} \left(\mathbf{t}_j^{\text{unit}} - (\hat{r}_{ij}\cdot\mathbf{t}_j^{\text{unit}})\,\hat{r}_{ij}\right)$$

$$\frac{\partial V_{\text{barrier}}}{\partial r_{ij}} = -\frac{12\, V_{\text{barrier}}}{r_{ij}}$$

By Newton's third law, $\mathbf{F}_j = -\mathbf{F}_i$.

> **Note:** Because backbone neighbor forces ($i\pm1$, $j\pm1$) are neglected, energy is not exactly conserved in NVE. Use NVT or Langevin dynamics.

## Parameters

| Parameter | Role |
|-----------|------|
| `epsilon` | Well depth |
| `rhb` | Equilibrium HB distance |
| `sigma` | Radial and angular width |
| `A` | Repulsive barrier prefactor |
