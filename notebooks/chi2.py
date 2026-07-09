"""
chi2.py — Cosmology Summer School
===================================
Likelihood functions for cosmological observables.

Students are encouraged to explore and modify this file
when adding new observables or changing the analysis.

Available functions:
    chi2_sne(Om, rho_de, z, mb, dmb)            — SNe Ia marginalised chi2
    chi2_sne_full(Om, Mcurl, rho_de, z, mb, dmb) — SNe Ia full chi2
    chi2_cc(Om, h, rho_de, z, H_obs, dH_obs)    — cosmic chronometers chi2
    chi2_bao(Om, h, rho_de, z_th, d_c_th, ...)  — DESI DR2 BAO chi2
"""

import numpy as np
from utils import solve_distances, DH, E2, z_to_a


# ---------------------------------------------------------------------------
# Type Ia Supernovae
# ---------------------------------------------------------------------------

def load_sne(dataset='pantheon+', data_dir='../data'):
    """Load SNe Ia data.

    Parameters
    ----------
    dataset  : str — 'pantheon+' (default), 'pantheon', or 'binned'
    data_dir : str — path to data directory

    Returns
    -------
    z   : array — redshifts (zHD for pantheon+, zcmb for others)
    mb  : array — corrected apparent magnitudes
    dmb : array — diagonal magnitude uncertainties
    """
    if dataset == 'pantheon+':
        data = np.genfromtxt(f'{data_dir}/Pantheon+SH0ES.dat',
                             names=True, dtype=None, encoding=None)
        # Hubble flow only: drop calibrators and very low-z
        mask = (data['IS_CALIBRATOR'] == 0) & (data['zHD'] > 0.01)
        z   = data['zHD'][mask]
        mb  = data['m_b_corr'][mask]
        dmb = data['m_b_corr_err_DIAG'][mask]

    elif dataset in ('pantheon', 'binned'):
        fname = 'pantheon.txt' if dataset == 'pantheon' else 'binned.txt'
        z, mb, dmb = np.loadtxt(f'{data_dir}/{fname}',
                                 usecols=(1, 4, 5), unpack=True, skiprows=1)
        idx = np.argsort(z)
        z, mb, dmb = z[idx], mb[idx], dmb[idx]

    else:
        raise ValueError(f"Unknown dataset '{dataset}'. Choose 'pantheon+', 'pantheon', or 'binned'.")

    return z, mb, dmb

def chi2_sne(Om, rho_de, z, mb, dmb):
    """Marginalised chi^2 for Type Ia supernovae.

    Analytically marginalises over script-M = M + 25 + 5*log10(c/H0),
    which absorbs both the absolute magnitude M and h completely.
    Only Omega_m (through d_L) remains as a free parameter.

        chi^2 = A - B^2/C

    where Delta_i = mb_i - 5*log10(d_L(z_i)) and:
        A = sum( Delta_i^2 / sigma_i^2 )
        B = sum( Delta_i   / sigma_i^2 )
        C = sum( 1         / sigma_i^2 )

    Parameters
    ----------
    Om     : float — Omega_m
    rho_de : function of a — Omega_DE(a)
    z      : array — redshifts of supernovae
    mb     : array — observed apparent magnitudes
    dmb    : array — uncertainties on mb

    Returns
    -------
    chi2   : float — marginalised chi^2
    Mcurl  : float — best-fit script-M
    """
    z_th, _, d_L_th = solve_distances(Om, rho_de, z_max=z.max() * 1.01,
                                       n_points=1000)
    d_L_data = np.interp(z, z_th, d_L_th)
    mu_th    = 5 * np.log10(d_L_data)

    delta  = mb - mu_th
    A      = np.sum((delta / dmb)**2)
    B      = np.sum( delta / dmb**2)
    C      = np.sum(    1  / dmb**2)
    Mcurl  = B / C
    chi2   = A - B**2 / C

    return chi2, Mcurl


def chi2_sne_full(Om, Mcurl, rho_de, z, mb, dmb):
    """Full chi^2 for Type Ia supernovae with explicit script-M.

        chi^2 = sum( (mb - 5*log10(d_L) - Mcurl)^2 / sigma^2 )

    For grid scans: precompute d_L outside the Mcurl loop for speed.

    Parameters
    ----------
    Om     : float — Omega_m
    Mcurl  : float — script-M = M + 25 + 5*log10(c/H0)
    rho_de : function of a — Omega_DE(a)
    z      : array — redshifts of supernovae
    mb     : array — observed apparent magnitudes
    dmb    : array — uncertainties on mb

    Returns
    -------
    chi2   : float
    """
    z_th, _, d_L_th = solve_distances(Om, rho_de, z_max=z.max() * 1.01,
                                       n_points=1000)
    d_L_data = np.interp(z, z_th, d_L_th)
    mu_th    = 5 * np.log10(d_L_data) + Mcurl
    return np.sum(((mb - mu_th) / dmb)**2)


# ---------------------------------------------------------------------------
# Cosmic chronometers
# ---------------------------------------------------------------------------

def chi2_cc(Om, h, rho_de, z, H_obs, dH_obs):
    """Chi^2 for cosmic chronometers.

    H_th(z) = 100 * h * E(z)

    Parameters
    ----------
    Om      : float — Omega_m
    h       : float — dimensionless Hubble parameter
    rho_de  : function of a — Omega_DE(a)
    z       : array — redshifts
    H_obs   : array — observed H(z) in km/s/Mpc
    dH_obs  : array — uncertainties on H(z)

    Returns
    -------
    chi2    : float
    """
    a    = z_to_a(z)
    H_th = 100 * h * np.sqrt(E2(a, Om, rho_de))
    return np.sum(((H_th - H_obs) / dH_obs)**2)


# ---------------------------------------------------------------------------
# BAO — DESI DR2
# ---------------------------------------------------------------------------

def chi2_bao(Om, h, rho_de, z_th, d_c_th,
             z_DM, DM_obs, dDM,
             z_DH, DH_obs, dDH,
             z_DV, DV_obs, dDV):
    """Chi^2 for BAO observables (DESI DR2).

    Distances in Mpc, divided by rd from the fitting formula:
        rd = 147.05 * (Om * h^2 / 0.1432)^(-0.32)  [Mpc]

    Observables:
        DM/rd — comoving distance / sound horizon
        DH/rd — Hubble distance c/H(z) / sound horizon
        DV/rd — volume distance [z*DM^2*DH]^(1/3) / sound horizon

    Parameters
    ----------
    Om, h    : float — cosmological parameters
    rho_de   : function of a — Omega_DE(a)
    z_th     : array — redshift grid from solve_distances
    d_c_th   : array — comoving distance (units of c/H0) on z_th grid
    z_DM     : array — redshifts of DM/rd measurements
    DM_obs   : array — observed DM/rd values
    dDM      : array — uncertainties on DM/rd
    z_DH     : array — redshifts of DH/rd measurements
    DH_obs   : array — observed DH/rd values
    dDH      : array — uncertainties on DH/rd
    z_DV     : array — redshifts of DV/rd measurements
    DV_obs   : array — observed DV/rd values
    dDV      : array — uncertainties on DV/rd

    Returns
    -------
    chi2     : float
    """
    # sound horizon [Mpc]
    rd = 147.05 * (Om * h**2 / 0.1432)**(-0.32)

    # Hubble distance c/H0 [Mpc]
    DH_h = DH(h)

    # interpolate comoving distance at data redshifts
    d_c_DM = np.interp(z_DM, z_th, d_c_th)
    d_c_DH = np.interp(z_DH, z_th, d_c_th)
    d_c_DV = np.interp(z_DV, z_th, d_c_th)

    # theory distances [Mpc]
    DM_th = d_c_DM * DH_h
    DH_th = DH_h / np.sqrt(E2(z_to_a(z_DH), Om, rho_de))
    DV_DM = d_c_DV * DH_h
    DV_DH = DH_h / np.sqrt(E2(z_to_a(z_DV), Om, rho_de))
    DV_th = (z_DV * DV_DM**2 * DV_DH)**(1/3)

    # chi2 contributions
    chi2  = np.sum(((DM_th/rd - DM_obs) / dDM)**2)
    chi2 += np.sum(((DH_th/rd - DH_obs) / dDH)**2)
    chi2 += np.sum(((DV_th/rd - DV_obs) / dDV)**2)

    return chi2
