"""
utils.py — Cosmology Summer School
===================================
Numerical routines for cosmological dynamics and distances.

Flat universe by default (Ok=0). Pass Ok != 0 to include curvature.
Radiation off by default (Or=0). Pass Or != 0 to include radiation.
Time is in units of H0^-1, distances in units of c/H0.
Use DH() and tH() to convert to physical units.

Radiation density today: Or = 4.18e-5 / h^2
    h=0.70 -> Or ~ 8.52e-5
    h=0.67 -> Or ~ 9.30e-5

Dark energy models available via de_model():
    'LCDM'  — cosmological constant (w = -1)
    'CPL'   — Chevallier-Polarski-Linder (w0, wa)
    'lCDM'  — lambda-CDM quintessence, Fonseca et al. 2021 (lam)
"""

import numpy as np
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
# Coordinate conversions
# ---------------------------------------------------------------------------

def a_to_z(a):
    """Scale factor to redshift: z = 1/a - 1."""
    return 1/a - 1

def z_to_a(z):
    """Redshift to scale factor: a = 1/(1+z)."""
    return 1 / (1 + z)


# ---------------------------------------------------------------------------
# Physical units
# ---------------------------------------------------------------------------

def DH(h=0.7):
    """Hubble distance c/H0 in Mpc.  DH ~ 4286 Mpc for h=0.7."""
    return 3e3 / h

def tH(h=0.7):
    """Hubble time 1/H0 in Gyr.  tH ~ 13.97 Gyr for h=0.7."""
    return 9.778 / h

def Or_obs(h=0.7):
    """Radiation density parameter today: Or = 4.18e-5 / h^2."""
    return 4.18e-5 / h**2


# ---------------------------------------------------------------------------
# Dark energy models
# ---------------------------------------------------------------------------

def de_model(model, Om, Ok=0, Or=0, **params):
    """Return (rho_de, drho_de) for a given dark energy model.

    OL = 1 - Om - Ok - Or  (derived from closure relation).

    Parameters
    ----------
    model  : str   — 'LCDM', 'CPL', or 'lCDM'
    Om     : float — Omega_m today
    Ok     : float — curvature parameter (default 0, flat)
    Or     : float — radiation parameter (default 0, no radiation)
                     use Or_obs(h) for the physical value
    **params : model-specific parameters
        CPL  : w0, wa
        lCDM : lam

    Returns
    -------
    rho_de  : function of a — Omega_DE(a)
    drho_de : function of a — d(Omega_DE)/dN,  N = ln(a)
    """
    OL = 1 - Om - Ok - Or

    if model == 'LCDM':
        rho_de  = lambda a: OL
        drho_de = lambda a: 0.0

    elif model == 'CPL':
        w0 = params['w0']
        wa = params['wa']
        rho_de  = lambda a: OL * a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
        drho_de = lambda a: rho_de(a) * (-3*(1+w0+wa) + 3*wa*a)

    elif model == 'lCDM':
        lam     = params['lam']
        rho_de  = lambda a: OL * a**(-lam**2)
        drho_de = lambda a: -lam**2 * OL * a**(-lam**2)

    else:
        raise ValueError(f"Unknown model '{model}'. Choose 'LCDM', 'CPL', or 'lCDM'.")

    return rho_de, drho_de


# ---------------------------------------------------------------------------
# Friedmann equation
# ---------------------------------------------------------------------------

def E2(a, Om, rho_de, Ok=0, Or=0):
    """Dimensionless Hubble parameter squared: E^2(a) = H^2(a)/H0^2.

        E^2(a) = Om/a^3 + Or/a^4 + Omega_DE(a) + Ok/a^2

    Radiation today: Or = 4.18e-5 / h^2  (e.g. ~8.5e-5 for h=0.7)
    """
    return Om / a**3 + Or / a**4 + rho_de(a) + Ok / a**2


def dE2_dN(a, Om, drho_de, Ok=0, Or=0):
    """Derivative of E^2 with respect to N = ln(a).

        dE^2/dN = -3*Om/a^3 - 4*Or/a^4 + d(Omega_DE)/dN - 2*Ok/a^2
    """
    return -3 * Om / a**3 - 4 * Or / a**4 + drho_de(a) - 2 * Ok / a**2


def adotdot(a, Om, rho_de, drho_de, Ok=0, Or=0):
    """Acceleration d^2a/dt^2 from the Raychaudhuri equation.

        adotdot = a * [ E2(a) + dE2_dN(a)/2 ]
    """
    return a * (E2(a, Om, rho_de, Ok, Or) + dE2_dN(a, Om, drho_de, Ok, Or) / 2)


# ---------------------------------------------------------------------------
# ODE solvers
# ---------------------------------------------------------------------------

def solve_friedmann(Om, rho_de, drho_de, Ok=0, Or=0,
                    t_forward=3.0, t_backward=2.0, n_points=2000):
    """Solve the Friedmann + Raychaudhuri system from today (a=1).

    State vector: y = [a, adot]

    Equations (dimensionless, H0=1):
        da/dt    = adot
        dadot/dt = a * [ E2(a) + dE2_dN(a)/2 ]

    Initial conditions:
        a(0)    = 1   (today)
        adot(0) = 1   (E(a=1) = 1 by definition)

    Parameters
    ----------
    Om         : float — Omega_m
    rho_de     : function of a — Omega_DE(a)
    drho_de    : function of a — d(Omega_DE)/dN
    Ok         : float — curvature parameter (default 0, flat)
    Or         : float — radiation parameter (default 0, no radiation)
    t_forward  : float — integration time into the future (H0^-1)
    t_backward : float — integration time into the past   (H0^-1)
    n_points   : int   — number of output points

    Returns
    -------
    t    : array — time in units of H0^-1
    a    : array — scale factor
    adot : array — da/dt
    """
    def rhs(t, y):
        a, adot = y
        if a <= 0:
            return [0.0, 0.0]
        add = adotdot(a, Om, rho_de, drho_de, Ok, Or)
        return [adot, add]

    # adot(today) = E(a=1) = 1 by definition of the Omega parameters
    y0 = [1.0, 1.0]

    def hit_zero(t, y):
        return y[0] - 1e-4
    hit_zero.terminal  = True
    hit_zero.direction = -1

    t_fwd = np.linspace(0,  t_forward,  n_points)
    t_bwd = np.linspace(0, -t_backward, n_points)

    sol_fwd = solve_ivp(rhs, [0,  t_forward],  y0, t_eval=t_fwd,
                        events=hit_zero, max_step=0.001, rtol=1e-8)
    sol_bwd = solve_ivp(rhs, [0, -t_backward], y0, t_eval=t_bwd,
                        events=hit_zero, max_step=0.001, rtol=1e-8)

    t    = np.concatenate([sol_bwd.t[::-1],    sol_fwd.t[1:]])
    a    = np.concatenate([sol_bwd.y[0][::-1], sol_fwd.y[0][1:]])
    adot = np.concatenate([sol_bwd.y[1][::-1], sol_fwd.y[1][1:]])

    mask = a > 0
    return t[mask], a[mask], adot[mask]


def solve_distances(Om, rho_de, Ok=0, Or=0, z_max=3.0, n_points=1000):
    """Compute comoving and luminosity distances by solving an ODE in z.

    Solves:
        d(d_c)/dz = 1 / E(z)

    with d_c(0) = 0. Returns the full distance vector at all redshifts
    from 0 to z_max in a single integration pass.

    Parameters
    ----------
    Om       : float — Omega_m
    rho_de   : function of a — Omega_DE(a)
    Ok       : float — curvature parameter (default 0, flat)
    Or       : float — radiation parameter (default 0, no radiation)
    z_max    : float — maximum redshift
    n_points : int   — number of output points

    Returns
    -------
    z   : array — redshift
    d_c : array — comoving distance  (units of c/H0)
    d_L : array — luminosity distance (units of c/H0)

    Multiply by DH(h) to convert to Mpc.
    """
    def Ez(z):
        return np.sqrt(np.abs(E2(z_to_a(z), Om, rho_de, Ok, Or)))

    def rhs(z, y):
        return [1.0 / Ez(z)]

    z_arr = np.linspace(0, z_max, n_points)

    sol = solve_ivp(rhs, [0, z_max], [0.0],
                    t_eval=z_arr, max_step=z_max/n_points, rtol=1e-8)

    z   = sol.t
    d_c = sol.y[0]
    d_L = (1 + z) * d_c

    return z, d_c, d_L
