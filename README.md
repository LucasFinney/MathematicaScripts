# Mathematica Scripts

A collection of Wolfram Mathematica scripts.

---

## HydrogenPlot.m

Generates 2D cross-sectional probability density plots for a single electron in a hydrogen-like atom ($Z = 1$).

### Physics

The electron wavefunction in spherical coordinates is:

$$\psi_{n,l,m}(r,\theta,\phi) = R_{n,l}(r)\, Y_l^m(\theta,\phi)$$

where $R_{n,l}$ is the radial wavefunction (built from generalized Laguerre polynomials) and $Y_l^m$ is the spherical harmonic. All quantities use **atomic units** ($a_0 = \hbar = m_e = e = 1$).

The plotted quantity is $|\psi_{n,l,m}|^2$ sampled in the $xz$-plane. Because $|\psi|^2 \propto |Y_l^m|^2$ and $|e^{im\phi}| = 1$, the density is azimuthally symmetric and any meridional cross-section is equivalent.

### Quantum Numbers

| Symbol | Name | Valid range |
|--------|------|-------------|
| `n` | Principal | $n \geq 1$ |
| `l` | Orbital angular momentum | $0 \leq l \leq n-1$ |
| `m` | Magnetic | $-l \leq m \leq l$ |
| `s` | Spin projection $m_s$ | $\pm 1/2$ |

Spin does not affect the spatial density; it is included in the plot label only.

### Usage

```mathematica
<< "HydrogenPlot.m"

(* Basic plot — 2p orbital, spin up *)
HydrogenPlot[2, 1, 0, 1/2]

(* Adjust gamma for brightness/contrast (default 0.3, range ~0.1–1.0) *)
HydrogenPlot[3, 2, 1, 1/2, "Gamma" -> 0.2]

(* Full quantum number table for n <= 4 *)
QuantumNumberTable[4]
```

The `"Gamma"` option applies a power-law transformation $|\psi|^{2\gamma}$ to the colour scale. Values below 1 compress the dynamic range, preventing the sharp nuclear peak from washing out orbital structure in the rest of the image.

### Output

Running the script produces:

1. A colour-coded table of all 60 valid $(n, l, m, m_s)$ states for $n \leq 4$, with Dirac ket notation.
2. A $4 \times 2$ grid of example density plots spanning the 1s through 4f orbitals.
