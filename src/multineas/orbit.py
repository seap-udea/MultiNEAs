"""
Orbital definitions and transformations for asteroids.

This module lets you define an asteroid orbit from:
- **Orbital elements** (q, e, i, Omega, w, M, a)
- **State vector** (position and velocity in Cartesian coordinates)
- **Impact parameters** (impact site, impact velocity, and date)

All transformation steps (state ↔ elements, geodetic ↔ ecliptic, impact → heliocentric state)
are exposed as standalone functions so you can use them independently. Units are documented
on each function: use consistent units (e.g. AU and AU/day with mu in AU³/day², or km and
km/s with mu in km³/s²). When in doubt, inputs and outputs use the same units.
"""

from __future__ import annotations

import os
import sys
from contextlib import contextmanager
from typing import Optional, Tuple, Union

import numpy as np
import rebound as rb
import spiceypy as spy
from astropy.time import Time
from matplotlib import pyplot as plt
from scipy.optimize import newton

from .plot import multineas_watermark
from .util import Util

__all__ = [
    # Constants
    "RAD_PER_DEG",
    "AU_KM",
    "DAY_S",
    "MU_SUN_AU",
    "MU_SUN_KM",
    # Unit conversion
    "state_km_per_s_to_au_per_day",
    "state_au_per_day_to_km_per_s",
    # State ↔ elements
    "transformation_x_to_e",
    "transformation_e_to_x",
    # Jacobian and period
    "compute_jacobian_x_to_e",
    "compute_orbital_period",
    # Geodetic / ecliptic
    "geo_to_rectangular",
    "geo_to_eclip",
    "get_velocity_ecliptic",
    "get_asteroid_state_vector",
    # Impact → orbit
    "get_orbit_from_impact",
    "compute_orbit_from_impact_integrated",
    "get_orbit_impactor",
    # Main class
    "Orbit",
    "OrbitalCoordinates",
    # Deprecated/Alias
    "get_pre_impact_orbital_elements",
    "get_state_vector",
]

# -----------------------------------------------------------------------------
# Compatibility class
# -----------------------------------------------------------------------------


class OrbitalCoordinates:
    """
    Compatibility class for coordinate transformations.

    This class wraps the standalone transformation functions in this module
    to provide backward compatibility with code that expects an OrbitalCoordinates
    instance (e.g. probability module).
    """

    def transformation_x_to_e(self, x, y, z, vx, vy, vz, mu):
        """Wrapper for transformation_x_to_e."""
        return transformation_x_to_e(x, y, z, vx, vy, vz, mu)

    def transformation_e_to_x(self, q, e, i, Omega, w, M, mu):
        """Wrapper for transformation_e_to_x."""
        return transformation_e_to_x(q, e, i, Omega, w, M, mu)

    def compute_jacobian_x_to_e(self, a, e, i, Omega, w, M, mu):
        """Wrapper for compute_jacobian_x_to_e."""
        return compute_jacobian_x_to_e(a, e, i, Omega, w, M, mu)

    def geo2eclip(self, lon, lat, alt, date=None, et=None):
        """
        Wrapper for geo_to_eclip (aliased as geo2eclip).
        """
        return geo_to_eclip(lon, lat, alt, date=date, et=et)

    def get_velocity_ecliptic(self, vx, vy, vz, lon, lat, alt, date=None, et=None):
        """Wrapper for get_velocity_ecliptic."""
        return get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date, et=et)


# -----------------------------------------------------------------------------
# Constants and unit conversion
# -----------------------------------------------------------------------------

# Angle: degrees to radians
RAD_PER_DEG: float = np.pi / 180.0

# Length: 1 AU in km (IAU 2012)
AU_KM: float = 149597870.7

# Time: 1 day in seconds
DAY_S: float = 86400.0

# Gravitational parameter of the Sun
# AU³/year² (Gaussian constant k², k = 0.01720209895)
MU_SUN_AU: float = 39
# km³/s²
MU_SUN_KM: float = 1.32712440018e11


def state_km_per_s_to_au_per_day(
    position_km: np.ndarray,
    velocity_km_s: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert state vector from km and km/s to AU and AU/day.

    Parameters
    ----------
    position_km : array-like, shape (3,)
        Position [km].
    velocity_km_s : array-like, shape (3,)
        Velocity [km/s].

    Returns
    -------
    position_au : np.ndarray, shape (3,)
        Position [AU].
    velocity_au_per_day : np.ndarray, shape (3,)
        Velocity [AU/day].
    """
    pos = np.asarray(position_km, dtype=float)
    vel = np.asarray(velocity_km_s, dtype=float)
    position_au = pos / AU_KM
    velocity_au_per_day = vel * DAY_S / AU_KM
    return position_au, velocity_au_per_day


def state_au_per_day_to_km_per_s(
    position_au: np.ndarray,
    velocity_au_per_day: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert state vector from AU and AU/day to km and km/s.

    Parameters
    ----------
    position_au : array-like, shape (3,)
        Position [AU].
    velocity_au_per_day : array-like, shape (3,)
        Velocity [AU/day].

    Returns
    -------
    position_km : np.ndarray, shape (3,)
        Position [km].
    velocity_km_s : np.ndarray, shape (3,)
        Velocity [km/s].
    """
    pos = np.asarray(position_au, dtype=float)
    vel = np.asarray(velocity_au_per_day, dtype=float)
    position_km = pos * AU_KM
    velocity_km_s = vel / DAY_S * AU_KM
    return position_km, velocity_km_s


@contextmanager
def _suppress_stdout_stderr() -> None:
    """
    Context manager to suppress stdout and stderr.

    This is used to hide REBOUND's verbose output when querying JPL Horizons.
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        try:
            sys.stdout = devnull
            sys.stderr = devnull
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


def transformation_x_to_e(
    x: float,
    y: float,
    z: float,
    vx: float,
    vy: float,
    vz: float,
    mu: float,
) -> Tuple[float, float, float, float, float, float, float]:
    """
    Transform state vector (Cartesian) to orbital elements.

    Converts position and velocity to (q, e, i, Omega, w, M, a) using SPICE.
    All inputs and outputs use the same unit system (see below).

    Parameters
    ----------
    x, y, z : float
        Position components [AU].
    vx, vy, vz : float
        Velocity components [AU/day].
    mu : float
        Gravitational parameter [AU³/day²].

    Returns
    -------
    q : float
        Periapsis distance [AU].
    e : float
        Eccentricity (dimensionless).
    i : float
        Inclination [radians].
    Omega : float
        Longitude of ascending node [radians].
    w : float
        Argument of periapsis [radians].
    M : float
        Mean anomaly [radians].
    a : float
        Semi-major axis [AU].

    Notes
    -----
    Units: position [AU], velocity [AU/day], mu [AU³/day²], angles [radians].
    Use `state_km_per_s_to_au_per_day` to convert from km and km/s first.
    """
    elements = spy.oscelt([x, y, z, vx, vy, vz], et=0, mu=mu)
    a = elements[0] / (1 - elements[1])
    return np.concatenate([elements[:6], [a]])


def transformation_e_to_x(
    q: float,
    e: float,
    i: float,
    Omega: float,
    w: float,
    M: float,
    mu: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Transform orbital elements to state vector (Cartesian).

    Inverse of `transformation_x_to_e`. Uses SPICE to compute position and
    velocity from (q, e, i, Omega, w, M).

    Parameters
    ----------
    q : float
        Periapsis distance [AU].
    e : float
        Eccentricity (dimensionless).
    i, Omega, w, M : float
        Inclination, longitude of node, argument of periapsis, mean anomaly
        [radians].
    mu : float
        Gravitational parameter [AU³/day²].

    Returns
    -------
    x, y, z : float
        Position [AU].
    vx, vy, vz : float
        Velocity [AU/day].

    Notes
    -----
    Units: same as inputs — lengths in AU, angles in radians, mu in AU³/day².
    """
    state_vec = spy.conics([q, e, i, Omega, w, M] + [0, mu], 0)

    r = state_vec[:3]
    v = state_vec[3:6]

    return r, v


def _kepler_equation(E, M, e):
    """
    Kepler's equation: M = E - e*sin(E)

    This function returns the residual: E - e*sin(E) - M = 0

    Parameters
    ----------
    E : float
        Eccentric anomaly (radians).
    M : float
        Mean anomaly (radians).
    e : float
        Eccentricity (dimensionless).

    Returns
    -------
    residual : float
        Residual of Kepler's equation: E - e*sin(E) - M
    """
    return E - e * np.sin(E) - M


def _sqrt_e(e: float) -> float:
    """
    Compute sqrt(1 - e^2).

    Parameters
    ----------
    e : float
        Eccentricity.

    Returns
    -------
    float
        sqrt(1 - e^2).
    """
    return (1 - e**2) ** 0.5


def _nu(a: float, mu: float) -> float:
    """
    Compute angular momentum per unit mass: sqrt(mu * a).

    Parameters
    ----------
    a : float
        Semi-major axis.
    mu : float
        Gravitational parameter.

    Returns
    -------
    float
        sqrt(mu * a).
    """
    return (mu * a) ** 0.5


def _r(a: float, e: float, E: float) -> float:
    """
    Compute radial distance from central body.

    Parameters
    ----------
    a : float
        Semi-major axis.
    e : float
        Eccentricity.
    E : float
        Eccentric anomaly.

    Returns
    -------
    float
        Radial distance: a * (1 - e * cos(E)).
    """
    return a * (1 - e * np.cos(E))


def _compute_functions(i: float, w: float, Omega: float, which: set) -> dict:
    """
    Compute transformation functions A, B, C, D, F, G.

    These are components of the rotation matrix that transforms from
    the orbital plane coordinate system to the inertial coordinate system.

    Parameters
    ----------
    i : float
        Inclination (radians).
    w : float
        Argument of periapsis (radians).
    Omega : float
        Longitude of ascending node (radians).
    which : set
        Set of strings indicating which functions to compute: {'A', 'B', 'C', 'D', 'F', 'G'}.

    Returns
    -------
    dict
        Dictionary containing the requested transformation functions.
    """
    results = {}

    if "A" in which:
        results["A"] = np.cos(Omega) * np.cos(w) - np.sin(Omega) * np.cos(i) * np.sin(w)
    if "B" in which:
        results["B"] = -np.cos(Omega) * np.sin(w) - np.sin(Omega) * np.cos(i) * np.cos(
            w
        )
    if "C" in which:
        results["C"] = np.sin(Omega) * np.cos(w) + np.cos(Omega) * np.cos(i) * np.sin(w)
    if "D" in which:
        results["D"] = -np.sin(Omega) * np.sin(w) + np.cos(Omega) * np.cos(i) * np.cos(
            w
        )
    if "F" in which:
        results["F"] = np.sin(w) * np.sin(i)
    if "G" in which:
        results["G"] = np.cos(w) * np.sin(i)

    return results


def _compute_state_vector(elements: list, mu: float) -> np.ndarray:
    """
    Compute state vector from orbital elements using SPICE.

    Parameters
    ----------
    elements : list
        Orbital elements [q, e, i, Omega, w, M].
    mu : float
        Gravitational parameter.

    Returns
    -------
    numpy.ndarray
        State vector [x, y, z, vx, vy, vz].
    """
    state_vector = spy.conics(elements + [0, mu], 0)
    return state_vector


def compute_jacobian_x_to_e(
    a: float, e: float, i: float, Omega: float, w: float, M: float, mu: float
) -> np.ndarray:
    """
    Compute the Jacobian matrix for transformation from Cartesian to orbital elements.

    Computes the Jacobian matrix J that relates differential changes in Cartesian
    coordinates (x, y, z, vx, vy, vz) to differential changes in orbital elements
    (q, e, i, Omega, w, M). The Jacobian is computed using analytical derivatives
    of the transformation equations.

    The Jacobian matrix has shape (6, 6) where:
    - Rows correspond to Cartesian coordinates: [x, y, z, vx, vy, vz]
    - Columns correspond to orbital elements: [q, e, i, Omega, w, M]

    Parameters
    ----------
    a : float
        Semi-major axis (AU).
    e : float
        Eccentricity (dimensionless).
    i : float
        Inclination (radians).
    Omega : float
        Longitude of ascending node (radians).
    w : float
        Argument of periapsis (radians).
    M : float
        Mean anomaly (radians).
    mu : float
        Standard gravitational parameter (AU^3/day^2).

    Returns
    -------
    J : numpy.ndarray
        Jacobian matrix of shape (6, 6). The matrix relates changes in Cartesian
        coordinates to changes in orbital elements (q, e, i, Omega, w, M).

    Examples
    --------
    >>> from multineas.orbit import OrbitalCoordinates
    >>> import numpy as np
    >>>
    >>> oc = OrbitalCoordinates()
    >>> # Define orbital elements
    >>> a = 1.0  # AU
    >>> e = 0.1
    >>> i = np.pi / 6  # 30 degrees
    >>> Omega = 0.0
    >>> w = 0.0
    >>> M = 0.0
    >>> mu = 0.01720209895**2
    >>>
    >>> # Compute Jacobian
    >>> J = oc.compute_jacobian_x_to_e(a, e, i, Omega, w, M, mu)
    >>> print(f"Jacobian shape: {J.shape}")
    >>> print(f"Jacobian determinant: {np.linalg.det(J):.6e}")

    Notes
    -----
    This method is particularly useful for computing probability density functions
    in orbital element space from probability densities in Cartesian space, as
    the Jacobian determinant gives the volume element transformation factor.

    The computation involves:
    1. Solving Kepler's equation to find the eccentric anomaly E
    2. Computing partial derivatives with respect to each orbital element
    3. Constructing the full Jacobian matrix
    4. Applying a transformation from (a, e) to (q, e) coordinates
    """
    # Solve Kepler's equation: M = E - e*sin(E) for E
    # Using M as initial guess (good for small eccentricities)
    E = newton(_kepler_equation, M, args=(M, e))

    functions = _compute_functions(i, w, Omega, {"A", "B", "C", "D", "F", "G"})
    A = functions["A"]
    B = functions["B"]
    C = functions["C"]
    D = functions["D"]
    F = functions["F"]
    G = functions["G"]

    r = _r(a, e, E)
    eps = _sqrt_e(e)
    nu = _nu(a, mu)
    nur = nu / r

    q = a * (1 - e)
    state_vector = _compute_state_vector([q, e, i, Omega, w, M], mu)

    # Partial derivative with respect to semi-major axis (a)
    partial_a_x = state_vector[0] / a
    partial_a_y = state_vector[1] / a
    partial_a_z = state_vector[2] / a
    partial_a_vx = -state_vector[3] / (2 * a)
    partial_a_vy = -state_vector[4] / (2 * a)
    partial_a_vz = -state_vector[5] / (2 * a)

    partial_a = [
        partial_a_x,
        partial_a_y,
        partial_a_z,
        partial_a_vx,
        partial_a_vy,
        partial_a_vz,
    ]

    # Partial derivative with respect to eccentricity (e)
    # dX/de
    dcosEde = -a * np.sin(E) ** 2 / r
    dsinEde = a * np.cos(E) * np.sin(E) / r
    dnurde = (nu * a / r**2) * (np.cos(E) - (a / r) * e * np.sin(E) ** 2)
    depsde = -e / eps

    drAde = a * (dcosEde - 1)
    drBde = a * (depsde * np.sin(E) + eps * dsinEde)

    dvAde = -(dnurde * np.sin(E) + nur * dsinEde)
    dvBde = dnurde * eps * np.cos(E) + nur * depsde * np.cos(E) + nur * eps * dcosEde

    partial_e = np.array(
        [
            drAde * A + drBde * B,
            drAde * C + drBde * D,
            drAde * F + drBde * G,
            dvAde * A + dvBde * B,
            dvAde * C + dvBde * D,
            dvAde * F + dvBde * G,
        ]
    )

    # Partial derivative with respect to inclination (i)
    partial_i_x = state_vector[2] * np.sin(Omega)
    partial_i_y = -state_vector[2] * np.cos(Omega)
    partial_i_z = -state_vector[0] * np.sin(Omega) + state_vector[1] * np.cos(Omega)

    partial_i_vx = state_vector[5] * np.sin(Omega)
    partial_i_vy = -state_vector[5] * np.cos(Omega)
    partial_i_vz = -state_vector[3] * np.sin(Omega) + state_vector[4] * np.cos(Omega)

    partial_i = [
        partial_i_x,
        partial_i_y,
        partial_i_z,
        partial_i_vx,
        partial_i_vy,
        partial_i_vz,
    ]

    # Partial derivative with respect to longitude of ascending node (Omega)
    partial_Omega_x = -state_vector[1]
    partial_Omega_y = state_vector[0]
    partial_Omega_z = 0

    partial_Omega_vx = -state_vector[4]
    partial_Omega_vy = state_vector[3]
    partial_Omega_vz = 0

    partial_Omega = [
        partial_Omega_x,
        partial_Omega_y,
        partial_Omega_z,
        partial_Omega_vx,
        partial_Omega_vy,
        partial_Omega_vz,
    ]

    # Partial derivative with respect to argument of periapsis (w)
    partial_w_x = -state_vector[1] * np.cos(i) - state_vector[2] * np.sin(i) * np.cos(
        Omega
    )
    partial_w_y = state_vector[0] * np.cos(i) - state_vector[2] * np.sin(i) * np.sin(
        Omega
    )
    partial_w_z = state_vector[0] * np.sin(i) * np.cos(Omega) + state_vector[
        1
    ] * np.sin(i) * np.sin(Omega)
    partial_w_vx = -state_vector[4] * np.cos(i) - state_vector[5] * np.sin(i) * np.cos(
        Omega
    )
    partial_w_vy = state_vector[3] * np.cos(i) - state_vector[5] * np.sin(i) * np.sin(
        Omega
    )
    partial_w_vz = state_vector[3] * np.sin(i) * np.cos(Omega) + state_vector[
        4
    ] * np.sin(i) * np.sin(Omega)

    partial_w = [
        partial_w_x,
        partial_w_y,
        partial_w_z,
        partial_w_vx,
        partial_w_vy,
        partial_w_vz,
    ]

    # Partial derivative with respect to mean anomaly (M)
    n = (mu / a**3) ** 0.5
    factor = -((mu * a**3) ** 0.5) / r**3

    partial_M_x = (1 / n) * state_vector[3]
    partial_M_y = (1 / n) * state_vector[4]
    partial_M_z = (1 / n) * state_vector[5]
    partial_M_vx = factor * state_vector[0]
    partial_M_vy = factor * state_vector[1]
    partial_M_vz = factor * state_vector[2]

    partial_M = [
        partial_M_x,
        partial_M_y,
        partial_M_z,
        partial_M_vx,
        partial_M_vy,
        partial_M_vz,
    ]

    # Construct the Jacobian matrix
    J = np.zeros((6, 6))
    J[:, 0] = partial_a
    J[:, 1] = partial_e
    J[:, 2] = partial_i
    J[:, 3] = partial_Omega
    J[:, 4] = partial_w
    J[:, 5] = partial_M

    # Transform from (a, e) to (q, e) coordinates
    Je2c = np.eye(6)
    Je2c[0, 0] = 1 / (1 - e)
    Je2c[0, 1] = q / (1 - e) ** 2
    JX2c = np.matmul(J, Je2c)

    return JX2c


def compute_orbital_period(a: float, mu: float) -> float:
    """
    Compute orbital period using Kepler's third law.

    Calculates the orbital period of a body in a Keplerian orbit using
    Kepler's third law: T = 2π * sqrt(a³/μ), where T is the period,
    a is the semi-major axis, and μ is the gravitational parameter.

    Parameters
    ----------
    a : float
        Semi-major axis [km] or [AU] (must match units of mu).
    mu : float
        Gravitational parameter [km³/s²] or [AU³/day²] (must match units of a).

    Returns
    -------
    float
        Orbital period [s] or [day] (matches time units implied by mu).

    Examples
    --------
    >>> from multineas.orbit import compute_orbital_period
    >>> import numpy as np
    >>>
    >>> # Earth's orbit (approximately)
    >>> a = 1.0  # AU
    >>> mu = 0.01720209895**2  # AU³/day²
    >>>
    >>> period = compute_orbital_period(a, mu)
    >>> print(f"Orbital period: {period:.2f} days")
    >>> print(f"Orbital period: {period/365.25:.2f} years")
    >>>
    >>> # Chelyabinsk-like orbit
    >>> a_km = 1.73 * 149597870  # km
    >>> mu_km = 1.32712440018e11  # km³/s² (Sun)
    >>> period_s = compute_orbital_period(a_km, mu_km)
    >>> print(f"Orbital period: {period_s/86400:.2f} days")

    Notes
    -----
    The units of the result depend on the units of the gravitational parameter:
    - If mu is in [km³/s²], the period is in [s]
    - If mu is in [AU³/day²], the period is in [day]
    """
    period = 2 * np.pi * np.sqrt(a**3 / mu)
    return period


def geo_to_rectangular(lon: float, lat: float, alt: float) -> np.ndarray:
    """
    Convert geodetic coordinates (longitude, latitude, altitude) to rectangular
    coordinates in the Earth-fixed frame (ITRF93).

    Uses SPICE to account for Earth's ellipsoidal shape (equatorial and polar
    radii from kernel data).

    Parameters
    ----------
    lon : float
        Geodetic longitude [degrees]. Positive eastward.
    lat : float
        Geodetic latitude [degrees]. Positive northward.
    alt : float
        Altitude above Earth's reference spheroid [km].

    Returns
    -------
    np.ndarray, shape (3,)
        Position vector [x, y, z] in Earth-fixed frame (ITRF93) [km].
    """
    lon_rad = lon * RAD_PER_DEG
    lat_rad = lat * RAD_PER_DEG

    n, props = spy.bodvrd("399", "RADII", 3)
    RE_spice = props[0]
    RP_spice = props[2]
    f_spice = (RE_spice - RP_spice) / RE_spice

    r_earth_fixed = spy.georec(lon_rad, lat_rad, alt, RE_spice, f_spice)
    return r_earth_fixed


def geo_to_eclip(
    lon: float,
    lat: float,
    alt: float,
    date: Optional[str] = None,
    et: Optional[float] = None,
    frame: str = "ITRF93",
) -> np.ndarray:
    """
    Convert geodetic coordinates to ecliptic J2000 coordinates.

    Converts geographic coordinates (latitude, longitude, altitude) of an impact
    event on Earth to ecliptic J2000 coordinates. This transformation accounts
    for Earth's rotation and applies the coordinate transformation from Earth-fixed
    to the inertial ecliptic frame.

    Parameters
    ----------
    lon : float
        Geodetic longitude [degrees].
    lat : float
        Geodetic latitude [degrees].
    alt : float
        Altitude above Earth's reference spheroid [km].
    date : str, optional
        UTC date and time in format 'YYYY-MM-DD HH:MM:SS'. Must be provided
        if et is not provided.
    et : float, optional
        Ephemeris time (ET) in seconds past J2000. Must be provided if date
        is not provided.
    frame : str, default 'ITRF93'
        Earth-fixed reference frame. Default is 'ITRF93'.

    Returns
    -------
    numpy.ndarray
        Position vector [x, y, z] in Ecliptic J2000 frame [km].

    Examples
    --------
    >>> from multineas.orbit import OrbitalCoordinates
    >>>
    >>> oc = OrbitalCoordinates()
    >>> # Chelyabinsk impact
    >>> lon = 61.1
    >>> lat = 54.8
    >>> alt = 30.0
    >>> date = '2013-02-15 03:20:33'
    >>>
    >>> r_eclip = oc.geo2eclip(lon, lat, alt, date=date)
    >>> print(f"Ecliptic coordinates: {r_eclip} km")
    >>> print(f"Magnitude: {np.linalg.norm(r_eclip):.2f} km")

    Notes
    -----
    The function first converts geodetic coordinates to Earth-centered Cartesian
    coordinates in the ITRF93 frame, then applies a transformation matrix to convert
    from ITRF93 (Earth-fixed) to ECLIPJ2000 (inertial, ecliptic-based) frame.
    This transformation is necessary for orbital calculations, ensuring the position
    is in an inertial reference frame.
    """
    r_earth_fixed = geo_to_rectangular(lon, lat, alt)

    # Convert ephemeris time if date is provided
    if date is not None:
        if et is not None:
            raise ValueError("Provide either 'date' or 'et', not both")
        et = spy.utc2et(date)
    elif et is None:
        raise ValueError("Either 'date' or 'et' must be provided")

    # Transform from Earth-fixed to Ecliptic J2000 frame
    M_ecl = spy.pxform(frame, "ECLIPJ2000", et)
    r_earth_ecl = spy.mxv(M_ecl, r_earth_fixed)

    return r_earth_ecl


def get_velocity_ecliptic(
    vx: float,
    vy: float,
    vz: float,
    lon: float,
    lat: float,
    alt: float,
    date: Optional[str] = None,
    et: Optional[float] = None,
) -> np.ndarray:
    """
    Convert velocity vector from Earth-fixed to ecliptic J2000 coordinates.

    This function takes a velocity vector in Earth-fixed coordinates and converts it
    to ecliptic J2000 coordinates, accounting for Earth's rotation. The observed
    velocity is corrected by adding the contribution from Earth's rotation.

    Parameters
    ----------
    vx, vy, vz : float
        Velocity components in Earth-fixed coordinates [km/s].
    lon, lat : float
        Geodetic longitude and latitude of the observation point [degrees].
    alt : float
        Altitude above Earth's reference spheroid [km].
    date : str, optional
        UTC date and time in format 'YYYY-MM-DD HH:MM:SS'. Must be provided
        if et is not provided.
    et : float, optional
        Ephemeris time (ET) in seconds past J2000. Must be provided if date
        is not provided.

    Returns
    -------
    numpy.ndarray
        Velocity vector [vx, vy, vz] in Ecliptic J2000 coordinates [km/s].

    Examples
    --------
    >>> from multineas.orbit import OrbitalCoordinates
    >>>
    >>> oc = OrbitalCoordinates()
    >>> # Chelyabinsk impact velocity
    >>> vx, vy, vz = 3.5, -12.8, -6.3  # km/s
    >>> lon, lat, alt = 61.1, 54.8, 30.0
    >>> date = '2013-02-15 03:20:33'
    >>>
    >>> v_eclip = oc.get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date)
    >>> print(f"Ecliptic velocity: {v_eclip} km/s")

    Notes
    -----
    The function accounts for Earth's rotation by adding the cross product
    of Earth's angular velocity and the position vector to the input velocity.
    This correction is necessary because the observed velocity is measured in
    the Earth-fixed frame, but orbital calculations require velocities in an
    inertial frame.
    """
    # Input validation
    if date is not None and et is not None:
        raise ValueError("Provide either 'date' or 'et', not both")
    if date is None and et is None:
        raise ValueError("Either 'date' or 'et' must be provided")

    # Convert velocity to numpy array
    v = np.array([vx, vy, vz])

    # Get position vector in Earth-fixed coordinates
    r = geo_to_rectangular(lon, lat, alt)

    # Earth's rotation parameters
    t_sidereal = 86164.09053083288  # Sidereal day in seconds
    w_earth = 2 * np.pi / t_sidereal  # Earth's angular velocity [rad/s]
    omega = np.array([0, 0, w_earth])

    # Add Earth's rotation contribution to velocity
    v_E = v + spy.vcrss(omega, r)  # Velocity in Earth-fixed frame [km/s]

    # Convert ephemeris time if date is provided
    if date is not None:
        et = spy.utc2et(date)

    # Transform from Earth-fixed to ecliptic J2000 coordinates
    mx = spy.pxform("ITRF93", "ECLIPJ2000", et)
    v_eclip = spy.mxv(mx, v_E)

    return v_eclip


def get_asteroid_state_vector(
    r_eclip: np.ndarray, v_eclip: np.ndarray, date: str
) -> np.ndarray:
    """
    Heliocentric state from geocentric impact position and velocity.

    Given the impact site position and velocity in Ecliptic J2000 (geocentric),
    adds Earth's position and velocity at the given date to obtain the asteroid's
    heliocentric state.

    Parameters
    ----------
    r_eclip : array-like, shape (3,)
        Geocentric position at impact in Ecliptic J2000 [km].
    v_eclip : array-like, shape (3,)
        Geocentric velocity at impact in Ecliptic J2000 [km/s].
    date : str
        UTC date/time, e.g. 'YYYY-MM-DD HH:MM:SS'.

    Returns
    -------
    np.ndarray, shape (6,)
        Heliocentric state (x, y, z, vx, vy, vz) in Ecliptic J2000 [km], [km/s].
    """
    rb.horizons.SSL_CONTEXT = "unverified"

    sim = rb.Simulation()
    sim.units = "km", "s", "kg"
    sim.integrator = "IAS15"
    sim.dt = -86400

    time = Time(date, format="iso")

    # Suppress REBOUND's verbose output when querying Horizons
    with _suppress_stdout_stderr():
        for i in ["Sun", "199", "299", "399", "499", "599", "699", "799", "899"]:
            sim.add(i, hash=f"{i}", date=f"JD{time.tdb.jd}")

    r_earth = np.array(sim.particles["399"].xyz)
    v_earth = np.array(sim.particles["399"].vxyz)

    r_asteroid = r_eclip + r_earth
    v_asteroid = v_eclip + v_earth  # así si es

    return np.concatenate([r_asteroid, v_asteroid])


def get_pre_impact_orbit(
    lon: float,
    lat: float,
    alt: float,
    vx: float,
    vy: float,
    vz: float,
    date: str,
    mu: float,
) -> Tuple[float, float, float, float, float, float, float]:
    """
    Compute pre-impact orbital elements at the moment of impact.

    Returns the *osculating* orbital elements at the time of impact: builds the
    heliocentric state from the impact site (geodetic), impact velocity
    (Earth-fixed), and date, then converts that state to elements. No N-body
    integration is performed. For the orbit obtained by integrating backward
    in time with REBOUND (Sun + planets), use
    `compute_orbit_from_impact_integrated` instead.

    All inputs are in physical units (km, km/s, degrees); output elements use
    AU and radians (mu must be in AU³/day²).

    Parameters
    ----------
    lon : float
        Geodetic longitude [degrees].
    lat : float
        Geodetic latitude [degrees].
    alt : float
        Altitude above Earth's reference spheroid [km].
    vx, vy, vz : float
        Velocity components in Earth-fixed coordinates [km/s].
    date : str
        UTC date and time, e.g. 'YYYY-MM-DD HH:MM:SS'.
    mu : float
        Gravitational parameter of the central body [AU³/day²]. Use `MU_SUN_AU`
        for the Sun.

    Returns
    -------
    q : float
        Periapsis distance [AU].
    e : float
        Eccentricity (dimensionless).
    i : float
        Inclination [radians].
    Omega : float
        Longitude of ascending node [radians].
    w : float
        Argument of periapsis [radians].
    M : float
        Mean anomaly [radians].
    a : float
        Semi-major axis [AU].
    """
    r_eclip = geo_to_eclip(lon, lat, alt, date=date)
    v_eclip = get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date)
    state_km = get_asteroid_state_vector(r_eclip, v_eclip, date)
    position_km = state_km[:3]
    velocity_km_s = state_km[3:6]
    position_au, velocity_au_per_day = state_km_per_s_to_au_per_day(
        position_km, velocity_km_s
    )
    elements = transformation_x_to_e(
        float(position_au[0]),
        float(position_au[1]),
        float(position_au[2]),
        float(velocity_au_per_day[0]),
        float(velocity_au_per_day[1]),
        float(velocity_au_per_day[2]),
        mu,
    )
    q = elements[0]
    e = elements[1]
    i = elements[2]
    Omega = np.mod(elements[3], 2 * np.pi)
    w = np.mod(elements[4], 2 * np.pi)
    M = elements[5]
    a = elements[6]
    return np.array([q, e, i, Omega, w, M, a], dtype=float)


def get_orbit_impactor(
    lon: float,
    lat: float,
    alt: float,
    vx: float,
    vy: float,
    vz: float,
    date: str,
    mu: float,
) -> Tuple[float, float, float, float, float, float, float]:
    """
    Compute the orbit of the impactor via REBOUND N-body backward integration.

    Unlike `get_orbit_from_impact` (which returns osculating elements *at* the
    moment of impact), this function adds the impactor to a REBOUND simulation
    with the Sun and major planets, integrates *backward* in time by one
    orbital period, then returns the orbital elements from the integrated
    state. That gives an orbit representative of the impactor's path before
    the encounter, with planetary perturbations included.

    Parameters
    ----------
    lon : float
        Geodetic longitude [degrees].
    lat : float
        Geodetic latitude [degrees].
    alt : float
        Altitude above Earth's reference spheroid [km].
    vx, vy, vz : float
        Velocity components in Earth-fixed coordinates [km/s].
    date : str
        UTC date and time, e.g. 'YYYY-MM-DD HH:MM:SS'.
    mu : float
        Gravitational parameter of the central body [AU³/day²]. Use `MU_SUN_AU`
        for the Sun.

    Returns
    -------
    q : float
        Periapsis distance [AU].
    e : float
        Eccentricity (dimensionless).
    i : float
        Inclination [radians].
    Omega : float
        Longitude of ascending node [radians].
    w : float
        Argument of periapsis [radians].
    M : float
        Mean anomaly [radians].
    a : float
        Semi-major axis [AU].
    """
    # Pre-impact heliocentric state at impact time (km, km/s)
    r_eclip = geo_to_eclip(lon, lat, alt, date=date)
    v_eclip = get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date)
    state_km = get_asteroid_state_vector(r_eclip, v_eclip, date)

    # Orbital period from osculating elements at impact (for integration time)
    position_au, velocity_au_per_day = state_km_per_s_to_au_per_day(
        state_km[:3], state_km[3:6]
    )
    elements = transformation_x_to_e(
        float(position_au[0]),
        float(position_au[1]),
        float(position_au[2]),
        float(velocity_au_per_day[0]),
        float(velocity_au_per_day[1]),
        float(velocity_au_per_day[2]),
        mu,
    )
    a_impact = elements[6]
    period_days = compute_orbital_period(a_impact, mu)
    period_s = period_days * DAY_S

    # REBOUND: Sun + planets at date, add impactor, integrate backward
    rb.horizons.SSL_CONTEXT = "unverified"
    sim = rb.Simulation()
    sim.units = "km", "s", "kg"
    sim.integrator = "IAS15"
    sim.dt = -period_s / 10.0

    time = Time(date, format="iso")
    with _suppress_stdout_stderr():
        for body in ["Sun", "199", "299", "399", "499", "599", "699", "799", "899"]:
            sim.add(body, hash=body, date=f"JD{time.tdb.jd}")

    sim.add(
        x=state_km[0],
        y=state_km[1],
        z=state_km[2],
        vx=state_km[3],
        vy=state_km[4],
        vz=state_km[5],
        hash="asteroid",
    )
    sim.integrate(-period_s)
    sim.move_to_hel()

    o = sim.particles["asteroid"].orbit()
    # REBOUND orbit: a in km, angles in radians; o.f = true anomaly
    a_au = o.a / AU_KM
    q = a_au * (1 - o.e)
    M = Util.true_anomaly_to_mean_anomaly(o.e, o.f)

    e = o.e
    i = o.inc
    Omega = np.mod(o.Omega, 2 * np.pi)
    w = np.mod(o.omega, 2 * np.pi)

    return np.array([q, e, i, Omega, w, M, a_au], dtype=float)


class Orbit:
    """
    Asteroid orbit defined from elements, state vector, or impact parameters.

    Internally the orbit is stored as orbital elements (q, e, i, Omega, w, M, a)
    and state vector (x, y, z, vx, vy, vz) in **AU** and **AU/day**, with angles
    in **radians**. The gravitational parameter `mu` must be in **AU³/day²**
    (e.g. use `MU_SUN_AU` for the Sun).

    You can construct an orbit in three ways:

    1. **Orbital elements**: pass `elements=(q, e, i, Omega, w, M, a)` in AU and radians.
    2. **State vector**: pass `state_vector=(x, y, z, vx, vy, vz)` in AU and AU/day.
    3. **Impact parameters**: pass `impact_parameters=(lon, lat, alt, vx, vy, vz, date)`
       with lon/lat in degrees, alt in km, velocity in km/s, date as UTC string.

    All transformation functions used under the hood are available as standalone
    functions in this module (e.g. `transformation_x_to_e`, `geo_to_eclip`) so you
    can reuse them with explicit control over units.
    """

    def __init__(
        self,
        mu: float,
        *,
        state_vector: Optional[Union[np.ndarray, Tuple[float, ...]]] = None,
        elements: Optional[Union[np.ndarray, Tuple[float, ...]]] = None,
        impact_parameters: Optional[Union[np.ndarray, Tuple]] = None,
    ) -> None:
        """
        Initialize an orbit from one of: state vector, elements, or impact parameters.

        Parameters
        ----------
        mu : float
            Gravitational parameter of the central body [AU³/day²]. Use
            `MU_SUN_AU` for heliocentric orbits.
        state_vector : array-like, optional
            (x, y, z, vx, vy, vz) in AU and AU/day. Exactly one of
            state_vector, elements, or impact_parameters must be given.
        elements : array-like, optional
            (q, e, i, Omega, w, M, a) in AU and radians. Exactly one of
            state_vector, elements, or impact_parameters must be given.
        impact_parameters : array-like, optional
            (lon, lat, alt, vx, vy, vz, date): lon/lat [deg], alt [km],
            velocity [km/s] in Earth-fixed frame, date as UTC string
            (e.g. 'YYYY-MM-DD HH:MM:SS'). Exactly one of state_vector,
            elements, or impact_parameters must be given.
        """
        self._mu = float(mu)
        self._impact_parameters: Optional[np.ndarray] = None

        if state_vector is not None:
            if elements is not None or impact_parameters is not None:
                raise ValueError(
                    "Provide exactly one of state_vector, elements, or impact_parameters."
                )
            vec = np.asarray(state_vector, dtype=float)
            if vec.size != 6:
                raise ValueError("state_vector must have 6 elements (x,y,z,vx,vy,vz).")
            self._state_vector = vec.ravel()
            self._elements = transformation_x_to_e(
                self._state_vector[0],
                self._state_vector[1],
                self._state_vector[2],
                self._state_vector[3],
                self._state_vector[4],
                self._state_vector[5],
                self._mu,
            )

        elif elements is not None:
            if impact_parameters is not None:
                raise ValueError(
                    "Provide exactly one of state_vector, elements, or impact_parameters."
                )
            el = np.asarray(elements, dtype=float)
            if el.size != 6:
                raise ValueError(
                    "elements must have 6 elements (q, e, i, Omega, w, M, a)."
                )
            self._elements = el.ravel()
            r, v = transformation_e_to_x(
                self._elements[0],
                self._elements[1],
                self._elements[2],
                self._elements[3],
                self._elements[4],
                self._elements[5],
                self._mu,
            )
            self._state_vector = np.concatenate([r, v])

        elif impact_parameters is not None:
            imp = impact_parameters
            if len(imp) < 7:
                raise ValueError(
                    "impact_parameters must have 7 elements: "
                    "lon, lat, alt, vx, vy, vz, date."
                )
            self._impact_parameters = np.array(imp, dtype=object)
            elements = get_orbit_impactor(
                float(imp[0]),
                float(imp[1]),
                float(imp[2]),
                float(imp[3]),
                float(imp[4]),
                float(imp[5]),
                str(imp[6]),
                self._mu,
            )
            self._elements = elements
            r, v = transformation_e_to_x(
                self._elements[0],
                self._elements[1],
                self._elements[2],
                self._elements[3],
                self._elements[4],
                self._elements[5],
                self._mu,
            )
            self._state_vector = np.concatenate([r, v])

        else:
            raise ValueError(
                "Provide exactly one of state_vector, elements, or impact_parameters."
            )

    @property
    def mu(self) -> float:
        """Gravitational parameter [AU³/day²]."""
        return self._mu

    @property
    def state_vector(self) -> np.ndarray:
        """State vector (x, y, z, vx, vy, vz) in AU and AU/day."""
        return self._state_vector.copy()

    @property
    def elements(self) -> np.ndarray:
        """Orbital elements (q, e, i, Omega, w, M, a) in AU and radians."""
        return self._elements.copy()

    @property
    def impact_parameters(self) -> Optional[np.ndarray]:
        """Impact parameters (lon, lat, alt, vx, vy, vz, date) or None."""
        return (
            None if self._impact_parameters is None else self._impact_parameters.copy()
        )

    @property
    def orbital_period(self) -> float:
        """Orbital period [day] from Kepler's third law."""
        if len(self._elements) == 6:
            a = self._elements[0] / (1 - self._elements[1])
        else:
            a = self._elements[6]
        return compute_orbital_period(a, self._mu)

    @property
    def jacobian_x_to_e(self) -> np.ndarray:
        """
        Jacobian of the transformation from state vector to orbital elements.

        Shape (6, 6): derivatives of (x,y,z,vx,vy,vz) with respect to (q,e,i,Omega,w,M).
        Units: same as state and elements (AU, AU/day, rad).
        """
        if len(self._elements) == 6:
            a = self._elements[0] / (1 - self._elements[1])
        else:
            a = self._elements[6]
        return compute_jacobian_x_to_e(
            a,
            self._elements[1],
            self._elements[2],
            self._elements[3],
            self._elements[4],
            self._elements[5],
            self._mu,
        )

    @property
    def pre_impact_orbit(self) -> Optional[np.ndarray]:
        """
        Orbital elements at impact (q, e, i, Omega, w, M, a) in AU and radians.

        None if this orbit was not defined from impact parameters.
        """
        if self._impact_parameters is None:
            return None
        elements = get_pre_impact_orbit(
            self._impact_parameters[0],
            self._impact_parameters[1],
            self._impact_parameters[2],
            self._impact_parameters[3],
            self._impact_parameters[4],
            self._impact_parameters[5],
            str(self._impact_parameters[6]),
        )
        return elements

    @property
    def pre_impact_state_vector(self) -> Optional[np.ndarray]:
        """
        State vector at impact (x, y, z, vx, vy, vz) in AU and AU/day.

        None if this orbit was not defined from impact parameters.
        """
        if self._impact_parameters is None:
            return None
        r_eclip = geo_to_eclip(
            self._impact_parameters[0],
            self._impact_parameters[1],
            self._impact_parameters[2],
            date=str(self._impact_parameters[6]),
        )
        v_eclip = get_velocity_ecliptic(
            self._impact_parameters[3],
            self._impact_parameters[4],
            self._impact_parameters[5],
            self._impact_parameters[0],
            self._impact_parameters[1],
            self._impact_parameters[2],
            date=str(self._impact_parameters[6]),
        )
        state_km = get_asteroid_state_vector(
            r_eclip, v_eclip, str(self._impact_parameters[6])
        )
        position_au, velocity_au_per_day = state_km_per_s_to_au_per_day(
            state_km[:3], state_km[3:6]
        )
        return np.concatenate([position_au, velocity_au_per_day])

    def plot(
        self,
        ax=None,
        *,
        n_points: int = 200,
        show_earth: bool = True,
        figsize: Tuple[float, float] = (6, 6),
    ):
        """
        Plot the orbit in the ecliptic plane (x–y) with Sun and optional Earth orbit.

        Similar to REBOUND orbit plots: Sun at origin, Earth on a circular orbit
        at 1 AU, and the asteroid orbit. All positions are in AU. Uses the
        MultiNEAs watermark and styles consistent with the plot module.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on. If None, a new figure and axes are created.
        n_points : int, optional
            Number of points to sample along the asteroid orbit. Default 200.
        show_earth : bool, optional
            If True, plot Earth's circular orbit (a = 1.0 AU). Default True.
        figsize : tuple, optional
            Figure size (width, height) in inches when ax is None. Default (8, 8).

        Returns
        -------
        fig : matplotlib.figure.Figure or None
            The figure (only if ax was None).
        ax : matplotlib.axes.Axes
            The axes used for the plot.

        Examples
        --------
        >>> from multineas.orbit2 import Orbit, MU_SUN_AU
        >>> import numpy as np
        >>> orb = Orbit(mu=MU_SUN_AU, elements=[0.7, 0.3, 0.1, 0, 0, 0, 1.0])
        >>> fig, ax = orb.plot()
        >>> ax.set_title("Heliocentric orbit")
        >>> plt.show()
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        else:
            fig = ax.get_figure()

        q, e, i, Omega, w, _M, _a = (
            self._elements[0],
            self._elements[1],
            self._elements[2],
            self._elements[3],
            self._elements[4],
            self._elements[5],
            self._elements[6],
        )

        # Sun at origin
        ax.scatter(
            [0],
            [0],
            s=180,
            c="#FDB813",
            edgecolors="#E5a00d",
            linewidths=1.5,
            zorder=5,
            label="Sun",
        )

        # Earth orbit (circular, a = 1.0 AU)
        if show_earth:
            theta_earth = np.linspace(0, 2 * np.pi, n_points)
            x_earth = np.cos(theta_earth)
            y_earth = np.sin(theta_earth)
            ax.plot(
                x_earth,
                y_earth,
                color="#2E86AB",
                lw=1.5,
                ls="-",
                label="Earth (a = 1.0 AU)",
                zorder=2,
            )

        # Asteroid orbit: sample by mean anomaly
        M_samples = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
        x_ast, y_ast = [], []
        for M_val in M_samples:
            r_i, _ = transformation_e_to_x(q, e, i, Omega, w, M_val, self._mu)
            x_ast.append(r_i[0])
            y_ast.append(r_i[1])
        x_ast = np.array(x_ast)
        y_ast = np.array(y_ast)
        ax.plot(
            x_ast,
            y_ast,
            color="#E94F37",
            lw=1.5,
            ls="-",
            label="Asteroid orbit",
            zorder=3,
        )

        # Current position on the orbit
        x0, y0 = self._state_vector[0], self._state_vector[1]
        ax.scatter(
            [x0],
            [y0],
            s=60,
            c="#E94F37",
            edgecolors="k",
            linewidths=0.8,
            marker="o",
            zorder=4,
            label="Asteroid (current)",
        )

        ax.axhline(0, color="k", lw=0.4, alpha=0.5, zorder=0)
        ax.axvline(0, color="k", lw=0.4, alpha=0.5, zorder=0)
        ax.set_aspect("equal")
        ax.set_xlabel("x [AU]", fontsize=11)
        ax.set_ylabel("y [AU]", fontsize=11)
        r_max = max(1.0, float(np.max(np.sqrt(x_ast**2 + y_ast**2))) * 1.15)
        ax.set_xlim(-r_max, r_max)
        ax.set_ylim(-r_max, r_max)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=9)
        multineas_watermark(ax, alpha=0.5)

        return fig, ax

    def __repr__(self) -> str:
        return (
            f"<Orbit q={self._elements[0]:.4f} AU e={self._elements[1]:.4f} "
            f"i={np.degrees(self._elements[2]):.2f}°>"
        )
    
# -----------------------------------------------------------------------------
# Aliases for backward compatibility
# -----------------------------------------------------------------------------
get_orbit_from_impact = get_pre_impact_orbit
get_pre_impact_orbital_elements = get_pre_impact_orbit
compute_orbit_from_impact_integrated = get_orbit_impactor
