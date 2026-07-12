# Quick Reference

Import everything with:
```python
from utils import *
from plots import *
from chi2 import chi2_sne, chi2_cc, chi2_bao, load_sne
```

---

## utils.py — numerics

### Coordinate conversions
```python
a_to_z(a)          # scale factor → redshift
z_to_a(z)          # redshift → scale factor
```

### Physical units
```python
DH(h=0.7)          # Hubble distance c/H0 in Mpc
tH(h=0.7)          # Hubble time 1/H0 in Gyr
Or_obs(h=0.7)      # radiation density today = 4.18e-5 / h^2
```

### Dark energy models
```python
rho_de, drho_de = de_model('LCDM',  Om)
rho_de, drho_de = de_model('CPL',   Om, w0=-0.9, wa=0.2)
rho_de, drho_de = de_model('lCDM',  Om, lam=0.5)

# with curvature or radiation (optional)
rho_de, drho_de = de_model('LCDM',  Om, Ok=-0.1)
rho_de, drho_de = de_model('LCDM',  Om, Or=Or_obs())
```

### Friedmann equation
```python
E2(a, Om, rho_de, Ok=0, Or=0)       # H^2(a) / H0^2
dE2_dN(a, Om, drho_de, Ok=0, Or=0)  # d(E^2)/d(ln a)
adotdot(a, Om, rho_de, drho_de, Ok=0, Or=0)  # d^2a/dt^2
```

### ODE solvers
```python
t, a, adot = solve_friedmann(Om, rho_de, drho_de,
                              Ok=0, Or=0,
                              t_forward=3.0, t_backward=2.0,
                              n_points=2000)

z, d_c, d_L = solve_distances(Om, rho_de,
                                Ok=0, Or=0,
                                z_max=3.0, n_points=1000)
```

---

## plots.py — visualisation

All plot functions accept optional range parameters.

```python
plot_at(t, a, Om, Ok=0, t_min=None, t_max=None)
plot_adot(t, adot, Om, Ok=0, t_min=None, t_max=None)
plot_phase(a, adot, Om, Ok=0, a_min=None, a_max=None)
plot_adotdot(t, a, Om, rho_de, drho_de, Ok=0, Or=0, t_min=None, t_max=None)
plot_densities(a, Om, rho_de, Ok=0, N_min=None, N_max=None)
plot_dL(z, d_L, Om, Ok=0, h=0.7, z_min=None, z_max=None)
plot_mu(z, d_L, Om, Ok=0, h=0.7, z_min=None, z_max=None)
```

---

## chi2.py — likelihoods

### Load data
```python
z, mb, dmb = load_sne('pantheon+')   # Pantheon+ (default)
z, mb, dmb = load_sne('pantheon')    # original Pantheon
z, mb, dmb = load_sne('binned')      # binned Pantheon (fast)
```

### Likelihoods
```python
# SNe Ia — marginalised over script-M
chi2, Mcurl = chi2_sne(Om, rho_de, z, mb, dmb)

# SNe Ia — full chi2 with explicit script-M
chi2 = chi2_sne_full(Om, Mcurl, rho_de, z, mb, dmb)

# cosmic chronometers
chi2 = chi2_cc(Om, h, rho_de, z, H_obs, dH_obs)

# BAO — pass precomputed distances for speed
z_th, d_c_th, _ = solve_distances(Om, rho_de, z_max=2.5)
chi2 = chi2_bao(Om, h, rho_de, z_th, d_c_th,
                z_DM, DM_obs, dDM,
                z_DH, DH_obs, dDH,
                z_DV, DV_obs, dDV)
```

---

## Dark energy models at a glance

| Model | `de_model` string | Extra params | ΛCDM limit |
|---|---|---|---|
| Cosmological constant | `'LCDM'` | — | — |
| CPL | `'CPL'` | `w0`, `wa` | `w0=-1, wa=0` |
| Fonseca et al. 2021 | `'lCDM'` | `lam` | `lam=0` |

Reference: Fonseca et al. 2021, arXiv:2104.14889

---

## Data files

| File | Contents |
|---|---|
| `Pantheon+SH0ES.dat` | 1701 Pantheon+ SNe Ia |
| `pantheon.txt` | 1048 original Pantheon SNe Ia |
| `cosmic_chrono.txt` | 32 cosmic chronometer $H(z)$ measurements |
| `desi_dr2_bao.txt` | 13 DESI DR2 BAO measurements |
