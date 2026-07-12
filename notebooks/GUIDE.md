# Cosmology Summer School — Notebook Guide

This guide walks through the four notebooks of the summer school, explaining the physics, the numerical choices, and how to explore each topic. Read it alongside the notebooks, not instead of them.

---

## Overview

The program follows a single arc:

1. **Dynamics** — how does the universe expand?
2. **Distances** — how do we convert expansion history into observable distances?
3. **Observables** — how do real datasets constrain the theory?
4. **Inference** — how do we extract parameter constraints from data?

Each notebook builds directly on the previous one. The numerical machinery lives in `utils.py`, plotting in `plots.py`, and likelihoods in `chi2.py` — you import these and focus on the physics.

---

## Notebook 1 — The expanding universe

### The Friedmann equation

The dynamics of the universe are governed by the Friedmann equation:

$$H^2(a) \equiv \left(\frac{\dot{a}}{a}\right)^2 = H_0^2\left[\frac{\Omega_m}{a^3} + \Omega_\Lambda + \frac{\Omega_k}{a^2}\right]$$

We work in dimensionless units throughout: time in units of $H_0^{-1}$ and distances in units of $c/H_0$. This means we never need to specify $H_0$ explicitly in the dynamics — it just sets the overall scale.

The key variable is the **scale factor** $a(t)$, normalised to $a = 1$ today. Redshift is related by $z = 1/a - 1$.

### The ODE system

Rather than solving the Friedmann equation as a single first-order ODE for $a$, we solve the **coupled system** $[a, \dot{a}]$:

$$\frac{da}{dt} = \dot{a}, \qquad \frac{d\dot{a}}{dt} = a\left[E^2(a) + \frac{1}{2}\frac{dE^2}{dN}\right]$$

where $N = \ln a$ and $E^2 = H^2/H_0^2$. This form, derived from the Friedmann and Raychaudhuri equations, handles recollapsing universes naturally — $\dot{a}$ simply changes sign at the turnaround without any intervention.

Initial conditions are $a(0) = 1$ and $\dot{a}(0) = 1$, which follow directly from the definition of the $\Omega$ parameters today.

### Things to explore

- **Curvature** — pass `Ok != 0` to `solve_friedmann` and `de_model`. Closed universes ($\Omega_k < 0$) can recollapse even with a positive cosmological constant if matter dominates enough.
- **Negative $\Lambda$** — set `OL < 0` in the ΛCDM model. The universe always recollapses.
- **The phase portrait** — `plot_phase(a, adot, Om)` shows $\dot{a}$ vs $a$. Zeros of $\dot{a}$ are turning points. The shape of this curve tells you the fate of the universe without integrating in time.
- **The acceleration** — `plot_adotdot` shows when the expansion transitions from deceleration to acceleration. For ΛCDM this happens at $z \approx 0.67$.
- **Energy densities** — `plot_densities` shows $\log|\rho_i|$ vs $N = \ln a$. The slopes are $-3$ for matter, $0$ for $\Lambda$, $-2$ for curvature — read off directly from the Friedmann equation.

---

## Notebook 2 — Cosmological distances

### From dynamics to distances

Once we have $H(z)$, all cosmological distances follow from a single integral — the **comoving distance**:

$$d_c(z) = \int_0^z \frac{dz'}{E(z')}$$

in units of $c/H_0$. Rather than evaluating this as a standard quadrature, we solve it as an ODE:

$$\frac{dd_c}{dz} = \frac{1}{E(z)}, \qquad d_c(0) = 0$$

This gives the full $d_c(z)$ array at all redshifts in a single pass — efficient and accurate.

The **luminosity distance** (relevant for SNe Ia) is:

$$d_L(z) = (1+z)\,d_c(z)$$

To convert to physical units (Mpc): multiply by $D_H(h) = c/H_0 = 3000/h$ Mpc.

### The distance modulus

Supernovae are characterised by their **distance modulus**:

$$\mu(z) = 5\log_{10}\left(d_L(z) \cdot D_H(h)\right) + 25$$

where $d_L$ is in units of $c/H_0$ and $D_H$ is in Mpc. The $+25$ converts from Mpc to the standard 10 pc reference distance.

### Dark energy models

Three models are available via `de_model(model, Om, **params)`:

| Model | String | Parameters | ΛCDM limit |
|---|---|---|---|
| Cosmological constant | `'LCDM'` | — | exact |
| Chevallier-Polarski-Linder | `'CPL'` | `w0`, `wa` | `w0=-1, wa=0` |
| Fonseca et al. 2021 | `'lCDM'` | `lam` | `lam=0` |

The Fonseca et al. model (arXiv:2104.14889) parametrises quintessence with a scalar field linear in $N = \ln a$, giving the analytic dark energy density $\Omega_{DE}(a) = \Omega_{DE,0}\,a^{-\lambda^2}$. No field evolution is needed — $H(z)$ is fully analytic.

### Things to explore

- Compare $d_L(z)$ for different models — the differences are small but measurable with SNe Ia out to $z \sim 2$.
- Try `lam = 0.5` in lCDM and compare to CPL with equivalent $w_0$, $w_a$.
- Notice that all models agree at $z = 0$ by construction ($d_L(0) = 0$) and diverge at high $z$.

---

## Notebook 3a — Type Ia Supernovae

### The observable

Type Ia supernovae are standard candles. Their observed apparent magnitude $m_b$ is related to theory by:

$$m_b = 5\log_{10}(d_L(z)) + \mathcal{M}$$

where $\mathcal{M} = M + 25 + 5\log_{10}(D_H)$ absorbs the absolute magnitude $M$ and the Hubble constant $H_0$. This means SNe alone cannot constrain $H_0$ — they constrain $\Omega_m$ and the combination $M + 5\log_{10}(H_0)$.

### Marginalisation

Since $\mathcal{M}$ is a nuisance parameter, we marginalise over it analytically. Writing $\Delta_i = m_b^i - 5\log_{10}(d_L(z_i))$, the marginalised $\chi^2$ is:

$$\chi^2 = A - \frac{B^2}{C}, \qquad A = \sum_i \frac{\Delta_i^2}{\sigma_i^2}, \quad B = \sum_i \frac{\Delta_i}{\sigma_i^2}, \quad C = \sum_i \frac{1}{\sigma_i^2}$$

The best-fit $\mathcal{M} = B/C$ is recovered automatically. From it, and an assumed absolute magnitude $M$, we can extract $H_0$:

$$H_0 = \frac{c}{10^{(\mathcal{M} - M - 25)/5}}$$

### The data

We use the **Pantheon+** sample (Brout et al. 2022, arXiv:2202.04077) — 1701 SNe Ia, of which 1624 are cosmological Hubble flow supernovae (excluding the 77 SH0ES Cepheid calibrators and applying a $z > 0.01$ cut). We use diagonal errors only, ignoring the full covariance matrix for simplicity.

---

## Notebook 3b — Cosmic Chronometers

### The observable

Cosmic chronometers measure $H(z)$ directly from the differential ages of passively evolving galaxies:

$$H(z) = -\frac{1}{1+z}\frac{dz}{dt}$$

The theoretical prediction is simply:

$$H_{th}(z) = 100\,h\,E(z)$$

Unlike SNe, $h$ enters directly — both $\Omega_m$ and $h$ are independently constrained. No distance integration is needed: `E2` is evaluated algebraically at the data redshifts.

### The data

We use the Moresco et al. compilation — 32 measurements from $z = 0.07$ to $z = 1.97$.

---

## Notebook 3c — Baryon Acoustic Oscillations

### The observable

The BAO scale — the sound horizon at the drag epoch $r_d \approx 147$ Mpc — appears as a characteristic scale in the galaxy distribution. By measuring it both transversely and radially, we access:

$$\frac{D_M(z)}{r_d}, \qquad \frac{D_H(z)}{r_d}, \qquad \frac{D_V(z)}{r_d}$$

where $D_M = d_c \cdot D_H$ is the comoving distance, $D_H(z) = c/H(z)$ is the Hubble distance, and $D_V = [z\,D_M^2\,D_H]^{1/3}$ is the angle-averaged distance used at low $z$ where the two cannot be separated.

### The sound horizon

We use the fitting formula (Eisenstein & Hu 1998):

$$r_d = 147.05\left(\frac{\Omega_m h^2}{0.1432}\right)^{-0.32} \text{ Mpc}$$

accurate to ~0.5%. The full computation requires integrating over the pre-recombination plasma, including radiation and baryons — a natural extension once the program is complete.

### The data

We use **DESI DR2** (Abdul Karim et al. 2025, arXiv:2503.14738) — 13 measurements from 7 tracer samples spanning $0.295 \leq z \leq 2.330$. Correlations between $D_M/r_d$ and $D_H/r_d$ at the same redshift are ignored for simplicity.

---

## Notebook 4 — MCMC Parameter Inference

### The posterior

The posterior distribution combines all datasets:

$$\ln\mathcal{P}(\theta) = \ln\pi(\theta) - \frac{1}{2}\sum_{\text{datasets}}\chi^2(\theta)$$

We use flat priors $\pi(\theta)$ within wide bounds. The datasets are switched on and off via simple boolean flags at the top of the notebook.

### The sampler

We use **emcee** (Foreman-Mackey et al. 2013) — an affine-invariant ensemble sampler that requires no covariance matrix tuning. Before running the MCMC, we strongly recommend working through the [emcee tutorial](https://emcee.readthedocs.io/en/stable/tutorials/line/) — it is one of the clearest introductions to ensemble MCMC available, using a simple line-fitting example that builds exactly the intuition needed here.

Key settings:
- **Walkers** — 32 is standard. More walkers help explore multimodal posteriors but are slower.
- **Steps** — 2000 is a starting point. Check the autocorrelation time $\tau$; you want at least $50\tau$ steps.
- **Burn-in** — discard the first 500 steps. Again, $\sim 2\tau$ is the guideline.

### Convergence

Check the **trace plots** — each walker should look like white noise after burn-in, with no drift or structure. The **autocorrelation time** $\tau$ is the number of steps between independent samples; the notebook computes it automatically.

### Corner plots

We use **getdist** for corner plots — the same tool used in published cosmological analyses including Planck, DESI, and Pantheon+. The triangle plot shows 1D marginalised posteriors on the diagonal and 2D joint posteriors off-diagonal, with $1\sigma$ and $2\sigma$ contours.

### Saving and loading chains

Chains are saved to `notebooks/chains/` as `.npy` files with accompanying `.json` metadata. Set `load_chain = True` in S7 and paste the chain filename to reproduce corner plots without rerunning the MCMC.

---

## Physical units

All internal calculations are dimensionless:

| Quantity | Unit |
|---|---|
| Time | $H_0^{-1}$ |
| Distance | $c/H_0$ |

To convert to physical units, use:

```python
DH(h)   # Hubble distance c/H0 in Mpc  (~4286 Mpc for h=0.7)
tH(h)   # Hubble time 1/H0 in Gyr      (~13.97 Gyr for h=0.7)
```

The radiation density today (off by default) is:

```python
Or_obs(h)   # = 4.18e-5 / h^2
```
