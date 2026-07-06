# Getting Started

Welcome to the Cosmology Summer School. This guide will get you up and running in a few steps.

---

## 1. Get the materials

Download or clone this repository:

```bash
git clone https://github.com/tgbar/2026-IA-summer.git
cd 2026-IA-summer
```

Or click the green **Code** button on GitHub and choose **Download ZIP**.

---

## 2. Install conda

We use conda to manage the Python environment. If you don't have it yet, install **Miniforge** (recommended):

- **Mac / Linux**: download and run the installer from https://github.com/conda-forge/miniforge
- **Windows**: same link, choose the Windows installer

Once installed, open a terminal and verify:

```bash
conda --version
```

---

## 3. Create the environment

From inside the repository folder:

```bash
mamba env create -f environment.yml
```

This installs Python and all the packages we need. It takes a few minutes the first time.

Then activate it:

```bash
conda activate cosmology
```

You should see `(cosmology)` at the start of your terminal prompt.

---

## 4. Open the notebooks

The simplest way is from the terminal:

```bash
jupyter lab
```

This opens JupyterLab in your browser. Navigate to the `notebooks/` folder and open `01_friedmann.ipynb` to start.

---

## 5. You are ready

Open `01_friedmann.ipynb`, select the `Python (cosmology)` kernel, and run the first cell.

---

## Optional: VS Code

If you prefer to work in VS Code instead of JupyterLab:

1. Download and install VS Code from https://code.visualstudio.com
2. Install the **Python** and **Jupyter** extensions (search in the Extensions panel)
3. Open the repository folder: `File → Open Folder`
4. Open a notebook and select `Python (cosmology)` as the kernel

---

## Optional: register the kernel for JupyterLab

If JupyterLab does not show the `cosmology` kernel, run this once:

```bash
python -m ipykernel install --user --name cosmology --display-name "Python (cosmology)"
```

Then restart JupyterLab.

---

## Packages included

| Package | Purpose |
|---|---|
| numpy | arrays and numerical computation |
| scipy | ODE solvers and integration |
| matplotlib | plotting |
| jupyterlab | notebook interface |
| ipywidgets | interactive widgets |
| emcee | MCMC sampling |
| getdist | posterior analysis and corner plots |
