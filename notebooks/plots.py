"""
plots.py — Cosmology Summer School
====================================
Plotting routines for cosmological dynamics and distances.
Imports numerical functions from utils.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from utils import adotdot, E2, dE2_dN

plt.rcParams.update({
    'figure.figsize': (8, 5),
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 2,
    'font.size': 12,
})


def _title(Om, Ok):
    """Standard title string showing Omega_m, Omega_Lambda, Omega_k."""
    OL = 1 - Om - Ok
    curvature = 'open' if Ok > 0 else 'closed' if Ok < 0 else 'flat'
    return f'$\\Omega_m={Om},\\ \\Omega_\\Lambda={OL:.2f},\\ \\Omega_k={Ok:.2f}$ — {curvature}'


# ---------------------------------------------------------------------------
# Dynamics plots
# ---------------------------------------------------------------------------

def plot_at(t, a, Om, Ok=0):
    """Plot scale factor a(t)."""
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, a, color='#2a78d6')
    ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.axhline(1, color='black', lw=0.5, alpha=0.3)
    ax.set_xlabel(r'time  $(H_0^{-1})$')
    ax.set_ylabel(r'scale factor  $a(t)$')
    ax.set_title(_title(Om, Ok))
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, a.max() * 1.1)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_adot(t, adot, Om, Ok=0):
    """Plot expansion rate adot(t)."""
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, adot, color='#2a78d6')
    ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.axhline(0, color='black', lw=0.8, alpha=0.3)
    ax.set_xlabel(r'time  $(H_0^{-1})$')
    ax.set_ylabel(r'$\dot{a}$')
    ax.set_title(_title(Om, Ok))
    ax.set_xlim(t[0], t[-1])
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_phase(a, adot, Om, Ok=0):
    """Plot phase portrait: adot vs a."""
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(a, adot, color='#e34948')
    ax.axhline(0, color='black', lw=0.8, alpha=0.3)
    ax.axvline(1, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.set_xlabel(r'scale factor  $a$')
    ax.set_ylabel(r'$\dot{a}$')
    ax.set_title(_title(Om, Ok))
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_adotdot(t, a, Om, rho_de, drho_de, Ok=0):
    """Plot acceleration adotdot(t)."""
    add = adotdot(a, Om, rho_de, drho_de, Ok)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, add, color='#4a3aa7')
    ax.axhline(0, color='black', lw=0.8, alpha=0.3)
    ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.set_xlabel(r'time  $(H_0^{-1})$')
    ax.set_ylabel(r'$\ddot{a}$')
    ax.set_title(_title(Om, Ok))
    ax.set_xlim(t[0], t[-1])
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_densities(a, Om, rho_de, Ok=0):
    """Plot log energy densities vs N = ln(a)."""
    OL = 1 - Om - Ok
    N_plot = np.linspace(np.log(a.min()), np.log(a.max()), 500)
    a_plot = np.exp(N_plot)

    rho_m = Om  * a_plot**(-3)
    rho_L = np.array([rho_de(ai) for ai in a_plot])
    rho_k = Ok  * a_plot**(-2)

    fig, ax = plt.subplots(figsize=(9, 5))

    def _plot_rho(rho, color, label):
        pos = rho >  0
        neg = rho <  0
        if np.any(pos):
            ax.plot(N_plot[pos], np.log10(rho[pos]),
                    color=color, ls='-',  label=label)
        if np.any(neg):
            ax.plot(N_plot[neg], np.log10(np.abs(rho[neg])),
                    color=color, ls=':', label=label + r' ($\rho<0$)')

    _plot_rho(rho_m, '#2a78d6', r'$\rho_m \propto a^{-3}$')
    _plot_rho(rho_L, '#e34948', r'$\rho_{DE}$')
    if Ok != 0:
        _plot_rho(rho_k, '#1baf7a', r'$\rho_k \propto a^{-2}$')

    ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.5, label='today ($N=0$)')
    ax.set_xlabel(r'$N = \ln a$')
    ax.set_ylabel(r'$\log_{10}|\rho_i|$')
    ax.set_title(_title(Om, Ok))
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Distance plots
# ---------------------------------------------------------------------------

def plot_dL(z, d_L, Om, Ok=0, h=0.7):
    """Plot luminosity distance d_L(z) in Mpc."""
    from utils import DH
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z, d_L * DH(h))
    ax.set_xlabel(r'redshift $z$')
    ax.set_ylabel(r'$d_L(z)$  [Mpc]')
    ax.set_title(_title(Om, Ok) + f', $h={h}$')
    plt.tight_layout()
    plt.show()


def plot_mu(z, d_L, Om, Ok=0, h=0.7):
    """Plot distance modulus mu(z)."""
    from utils import DH
    mu = 5 * np.log10(d_L * DH(h) * 1e6 / 10)  # d_L in pc
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z, mu)
    ax.set_xlabel(r'redshift $z$')
    ax.set_ylabel(r'distance modulus $\mu$  [mag]')
    ax.set_title(_title(Om, Ok) + f', $h={h}$')
    plt.tight_layout()
    plt.show()
