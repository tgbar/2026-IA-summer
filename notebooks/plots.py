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


def _t_mask(t, t_min, t_max):
    """Boolean mask for time range. None means use full range."""
    t_min = t_min if t_min is not None else t[0]
    t_max = t_max if t_max is not None else t[-1]
    return (t >= t_min) & (t <= t_max), t_min, t_max


def _a_mask(a, a_min, a_max):
    """Boolean mask for scale factor range. None means use full range."""
    a_min = a_min if a_min is not None else a[0]
    a_max = a_max if a_max is not None else a[-1]
    return (a >= a_min) & (a <= a_max), a_min, a_max


# ---------------------------------------------------------------------------
# Dynamics plots
# ---------------------------------------------------------------------------

def plot_at(t, a, Om, Ok=0, t_min=None, t_max=None):
    """Plot scale factor a(t).

    Parameters
    ----------
    t_min, t_max : float or None — time range to plot (H0^-1)
    """
    mask, t_min, t_max = _t_mask(t, t_min, t_max)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t[mask], a[mask], color='#2a78d6')
    ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.axhline(1, color='black', lw=0.5, alpha=0.3)
    ax.set_xlabel(r'time  $(H_0^{-1})$')
    ax.set_ylabel(r'scale factor  $a(t)$')
    ax.set_title(_title(Om, Ok))
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(0, a[mask].max() * 1.1)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_adot(t, adot, Om, Ok=0, t_min=None, t_max=None):
    """Plot expansion rate adot(t).

    Parameters
    ----------
    t_min, t_max : float or None — time range to plot (H0^-1)
    """
    mask, t_min, t_max = _t_mask(t, t_min, t_max)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t[mask], adot[mask], color='#2a78d6')
    ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.axhline(0, color='black', lw=0.8, alpha=0.3)
    ax.set_xlabel(r'time  $(H_0^{-1})$')
    ax.set_ylabel(r'$\dot{a}$')
    ax.set_title(_title(Om, Ok))
    ax.set_xlim(t_min, t_max)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_phase(a, adot, Om, Ok=0, a_min=None, a_max=None):
    """Plot phase portrait: adot vs a.

    Parameters
    ----------
    a_min, a_max : float or None — scale factor range to plot
    """
    mask, a_min, a_max = _a_mask(a, a_min, a_max)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(a[mask], adot[mask], color='#e34948')
    ax.axhline(0, color='black', lw=0.8, alpha=0.3)
    ax.axvline(1, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.set_xlabel(r'scale factor  $a$')
    ax.set_ylabel(r'$\dot{a}$')
    ax.set_title(_title(Om, Ok))
    ax.set_xlim(a_min, a_max)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_adotdot(t, a, Om, rho_de, drho_de, Ok=0, Or=0, t_min=None, t_max=None):
    """Plot acceleration adotdot(t).

    Parameters
    ----------
    t_min, t_max : float or None — time range to plot (H0^-1)
    """
    mask, t_min, t_max = _t_mask(t, t_min, t_max)
    add = adotdot(a[mask], Om, rho_de, drho_de, Ok, Or)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t[mask], add, color='#4a3aa7')
    ax.axhline(0, color='black', lw=0.8, alpha=0.3)
    ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.5, label='today')
    ax.set_xlabel(r'time  $(H_0^{-1})$')
    ax.set_ylabel(r'$\ddot{a}$')
    ax.set_title(_title(Om, Ok))
    ax.set_xlim(t_min, t_max)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


def plot_densities(a, Om, rho_de, Ok=0, a_min=None, a_max=None):
    """Plot log energy densities vs N = ln(a).

    Parameters
    ----------
    a_min, a_max : float or None — scale factor range to plot
    """
    _, a_min, a_max = _a_mask(a, a_min, a_max)
    N_plot = np.linspace(np.log(a_min), np.log(a_max), 500)
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

def plot_dL(z, d_L, Om, Ok=0, h=0.7, z_min=None, z_max=None):
    """Plot luminosity distance d_L(z) in Mpc.

    Parameters
    ----------
    z_min, z_max : float or None — redshift range to plot
    """
    from utils import DH
    z_min = z_min if z_min is not None else z[0]
    z_max = z_max if z_max is not None else z[-1]
    mask  = (z >= z_min) & (z <= z_max)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z[mask], d_L[mask] * DH(h))
    ax.set_xlabel(r'redshift $z$')
    ax.set_ylabel(r'$d_L(z)$  [Mpc]')
    ax.set_title(_title(Om, Ok) + f', $h={h}$')
    ax.set_xlim(z_min, z_max)
    plt.tight_layout()
    plt.show()


def plot_mu(z, d_L, Om, Ok=0, h=0.7, z_min=None, z_max=None):
    """Plot distance modulus mu(z).

    Parameters
    ----------
    z_min, z_max : float or None — redshift range to plot
    """
    from utils import DH
    z_min = z_min if z_min is not None else z[0]
    z_max = z_max if z_max is not None else z[-1]
    mask  = (z >= z_min) & (z <= z_max)
    mu = 5 * np.log10(d_L[mask] * DH(h)) + 25
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z[mask], mu)
    ax.set_xlabel(r'redshift $z$')
    ax.set_ylabel(r'distance modulus $\mu$  [mag]')
    ax.set_title(_title(Om, Ok) + f', $h={h}$')
    ax.set_xlim(z_min, z_max)
    plt.tight_layout()
    plt.show()
