# 2026 IA Cosmology Summer School

A three-week intensive program covering the dynamics of the universe, dark energy models, cosmological observables, and Bayesian parameter inference. 

---

## Getting started

See **`GETTING_STARTED.md`** for step by step installation instructions.

The conda environment with all required packages is in **`environment.yml`**.

For documentation on the notebooks see **`notebooks/GUIDE.md`**.
For a quick function reference see **`notebooks/REFERENCE.md`**.

---

## Notebooks

| File | Description |
|---|---|
| `01_friedmann.ipynb` | The expanding universe — scale factor, phase portrait, acceleration |
| `02_distances.ipynb` | Cosmological distances and the distance modulus |
| `03_sne.ipynb` | Type Ia supernovae — Pantheon+ data and chi2 analysis |
| `03_cc.ipynb` | Cosmic chronometers — H(z) measurements |
| `03_bao.ipynb` | Baryon acoustic oscillations — DESI DR2 data |
| `04_mcmc.ipynb` | MCMC parameter inference with emcee and getdist |

Supporting modules in `notebooks/`: `utils.py` (numerics), `plots.py` (plotting), `chi2.py` (likelihoods).

## Data

| File | Contents |
|---|---|
| `Pantheon+SH0ES.dat` | 1701 Pantheon+ Type Ia supernovae |
| `pantheon.txt` | 1048 original Pantheon supernovae |
| `cosmic_chrono.txt` | 32 cosmic chronometer H(z) measurements |
| `desi_dr2_bao.txt` | 13 DESI DR2 BAO measurements |

## Slides

Lecture slides in PDF covering the theory behind each notebook are in `slides/`.
