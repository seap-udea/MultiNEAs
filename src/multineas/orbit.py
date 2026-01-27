"""
Orbital coordinate transformation utilities.

This module provides tools for converting between Cartesian coordinates (position and velocity)
and orbital elements, as well as computing Jacobian matrices for these transformations.
It also includes coordinate transformations from geographic to ecliptic coordinates and
orbit integration capabilities using REBOUND.
"""

import numpy as np
from scipy.optimize import newton
import spiceypy as spy
from typing import Optional, Union, List, Tuple


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


class OrbitalCoordinates:
    """
    Class for transforming between Cartesian coordinates and orbital elements.
    
    This class provides methods to convert between position-velocity (X) coordinates
    and orbital elements (E), as well as computing the Jacobian matrix for these
    transformations. These transformations are essential for working with orbital
    mechanics and probability density functions in orbital element space.
    
    Attributes
    ----------
    None
    
    Examples
    --------
    Convert from Cartesian coordinates to orbital elements:
    
    >>> from multineas.orbit import OrbitalCoordinates
    >>> import numpy as np
    >>> 
    >>> # Define position and velocity in AU and AU/day
    >>> x, y, z = 1.0, 0.0, 0.0
    >>> vx, vy, vz = 0.0, 6.28, 0.0
    >>> mu = 0.01720209895**2  # Standard gravitational parameter
    >>> 
    >>> oc = OrbitalCoordinates()
    >>> q, e, i, Omega, w, M, a = oc.transformation_x_to_e(x, y, z, vx, vy, vz, mu)
    >>> print(f"Semi-major axis: {a:.4f} AU")
    >>> print(f"Eccentricity: {e:.4f}")
    
    Convert from orbital elements back to Cartesian coordinates:
    
    >>> x_new, y_new, z_new, vx_new, vy_new, vz_new = oc.transformation_e_to_x(
    ...     q, e, i, Omega, w, M, mu
    ... )
    >>> print(f"Position: ({x_new:.4f}, {y_new:.4f}, {z_new:.4f}) AU")
    
    Compute the Jacobian matrix for coordinate transformation:
    
    >>> J = oc.compute_jacobian_x_to_e(a, e, i, Omega, w, M, mu)
    >>> print(f"Jacobian shape: {J.shape}")
    """
    
    @staticmethod
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
        return (1 - e**2)**0.5
    
    @staticmethod
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
        return (mu * a)**0.5
    
    @staticmethod
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
    
    @staticmethod
    def _compute_functions(i: float, w: float, Omega: float, 
                          which: set) -> dict:
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
        
        if 'A' in which:
            results['A'] = (np.cos(Omega) * np.cos(w) - np.sin(Omega) * np.cos(i) * np.sin(w))
        if 'B' in which:
            results['B'] = (-np.cos(Omega) * np.sin(w) - np.sin(Omega) * np.cos(i) * np.cos(w))
        if 'C' in which:
            results['C'] = (np.sin(Omega) * np.cos(w) + np.cos(Omega) * np.cos(i) * np.sin(w))
        if 'D' in which:
            results['D'] = (-np.sin(Omega) * np.sin(w) + np.cos(Omega) * np.cos(i) * np.cos(w))
        if 'F' in which:
            results['F'] = np.sin(w) * np.sin(i)
        if 'G' in which:
            results['G'] = np.cos(w) * np.sin(i)
        
        return results
    
    @staticmethod
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
    
    def transformation_x_to_e(self, x: float, y: float, z: float, 
                             vx: float, vy: float, vz: float, 
                             mu: float) -> tuple[float, float, float, float, float, float, float]:
        """
        Transform from Cartesian coordinates (position and velocity) to orbital elements.
        
        Converts a state vector (position x, y, z and velocity vx, vy, vz) to orbital
        elements using SPICE routines. The orbital elements returned are: periapsis distance
        (q), eccentricity (e), inclination (i), longitude of ascending node (Omega),
        argument of periapsis (w), mean anomaly (M), and semi-major axis (a).
        
        Parameters
        ----------
        x : float
            X-component of position (AU).
        y : float
            Y-component of position (AU).
        z : float
            Z-component of position (AU).
        vx : float
            X-component of velocity (AU/day).
        vy : float
            Y-component of velocity (AU/day).
        vz : float
            Z-component of velocity (AU/day).
        mu : float
            Standard gravitational parameter (AU^3/day^2).
            
        Returns
        -------
        q : float
            Periapsis distance (AU).
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
        a : float
            Semi-major axis (AU).
            
        Examples
        --------
        >>> from multineas.orbit import OrbitalCoordinates
        >>> import numpy as np
        >>> 
        >>> oc = OrbitalCoordinates()
        >>> # Circular orbit at 1 AU
        >>> x, y, z = 1.0, 0.0, 0.0
        >>> vx, vy, vz = 0.0, 6.28, 0.0
        >>> mu = 0.01720209895**2
        >>> 
        >>> q, e, i, Omega, w, M, a = oc.transformation_x_to_e(x, y, z, vx, vy, vz, mu)
        >>> print(f"Eccentricity: {e:.6f}")  # Should be close to 0 for circular orbit
        >>> print(f"Semi-major axis: {a:.4f} AU")
        """
        elements = spy.oscelt([x, y, z, vx, vy, vz], et=0, mu=mu)
        q = elements[0]
        e = elements[1]
        i = elements[2]
        Omega = elements[3]
        w = elements[4]
        M = elements[5]
        a = q/(1-e)
        
        return q, e, i, Omega, w, M, a
    
    def transformation_e_to_x(self, q: float, e: float, i: float, 
                              Omega: float, w: float, M: float, 
                              mu: float) -> tuple[float, float, float, float, float, float]:
        """
        Transform from orbital elements to Cartesian coordinates (position and velocity).
        
        Converts orbital elements to a state vector (position x, y, z and velocity
        vx, vy, vz) using SPICE routines. This is the inverse transformation of
        `transformation_x_to_e`.
        
        Parameters
        ----------
        q : float
            Periapsis distance (AU).
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
        x : float
            X-component of position (AU).
        y : float
            Y-component of position (AU).
        z : float
            Z-component of position (AU).
        vx : float
            X-component of velocity (AU/day).
        vy : float
            Y-component of velocity (AU/day).
        vz : float
            Z-component of velocity (AU/day).
            
        Examples
        --------
        >>> from multineas.orbit import OrbitalCoordinates
        >>> import numpy as np
        >>> 
        >>> oc = OrbitalCoordinates()
        >>> # Define orbital elements
        >>> q = 0.5  # AU
        >>> e = 0.5
        >>> i = np.pi / 4  # 45 degrees
        >>> Omega = 0.0
        >>> w = 0.0
        >>> M = 0.0
        >>> mu = 0.01720209895**2
        >>> 
        >>> x, y, z, vx, vy, vz = oc.transformation_e_to_x(q, e, i, Omega, w, M, mu)
        >>> print(f"Position: ({x:.4f}, {y:.4f}, {z:.4f}) AU")
        >>> print(f"Velocity: ({vx:.4f}, {vy:.4f}, {vz:.4f}) AU/day")
        """
        state_vec = spy.conics([q, e, i, Omega, w, M]+[0, mu], 0)
        x = state_vec[0]
        y = state_vec[1]
        z = state_vec[2]
        vx = state_vec[3]
        vy = state_vec[4]
        vz = state_vec[5]
        
        return x, y, z, vx, vy, vz
    
    def compute_jacobian_x_to_e(self, a: float, e: float, i: float, 
                                Omega: float, w: float, M: float, 
                                mu: float) -> np.ndarray:
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
        
        functions = self._compute_functions(i, w, Omega, {'A', 'B', 'C', 'D', 'F', 'G'})
        A = functions['A']
        B = functions['B']
        C = functions['C']
        D = functions['D']
        F = functions['F']
        G = functions['G']
        
        r = self._r(a, e, E)
        eps = self._sqrt_e(e)
        nu = self._nu(a, mu)
        nur = nu/r
        
        q = a*(1-e)
        state_vector = self._compute_state_vector([q, e, i, Omega, w, M], mu)
        
        # Partial derivative with respect to semi-major axis (a)
        partial_a_x = state_vector[0]/a
        partial_a_y = state_vector[1]/a
        partial_a_z = state_vector[2]/a
        partial_a_vx = -state_vector[3]/(2*a)
        partial_a_vy = -state_vector[4]/(2*a)
        partial_a_vz = -state_vector[5]/(2*a)
        
        partial_a = [partial_a_x, partial_a_y, partial_a_z, partial_a_vx, partial_a_vy, partial_a_vz]
        
        # Partial derivative with respect to eccentricity (e)
        # dX/de
        dcosEde = -a*np.sin(E)**2/r       
        dsinEde = a*np.cos(E)*np.sin(E)/r
        dnurde = (nu*a/r**2)*(np.cos(E)-(a/r)*e*np.sin(E)**2)
        depsde = -e/eps
        
        drAde = a*(dcosEde-1)
        drBde = a*(depsde*np.sin(E)+eps*dsinEde)
        
        dvAde = -(dnurde*np.sin(E)+nur*dsinEde)
        dvBde = (dnurde*eps*np.cos(E)+nur*depsde*np.cos(E)+nur*eps*dcosEde)
        
        partial_e = np.array([
            drAde*A+drBde*B,
            drAde*C+drBde*D,
            drAde*F+drBde*G,
            dvAde*A+dvBde*B,
            dvAde*C+dvBde*D,
            dvAde*F+dvBde*G
        ])
        
        # Partial derivative with respect to inclination (i)
        partial_i_x = state_vector[2]*np.sin(Omega)
        partial_i_y = -state_vector[2]*np.cos(Omega)
        partial_i_z = -state_vector[0]*np.sin(Omega) + state_vector[1]*np.cos(Omega)
        
        partial_i_vx = state_vector[5]*np.sin(Omega)
        partial_i_vy = -state_vector[5]*np.cos(Omega)
        partial_i_vz = -state_vector[3]*np.sin(Omega) + state_vector[4]*np.cos(Omega)
        
        partial_i = [partial_i_x, partial_i_y, partial_i_z, partial_i_vx, partial_i_vy, partial_i_vz]
        
        # Partial derivative with respect to longitude of ascending node (Omega)
        partial_Omega_x = -state_vector[1]
        partial_Omega_y = state_vector[0]
        partial_Omega_z = 0
        
        partial_Omega_vx = -state_vector[4]
        partial_Omega_vy = state_vector[3]
        partial_Omega_vz = 0
        
        partial_Omega = [partial_Omega_x, partial_Omega_y, partial_Omega_z, 
                        partial_Omega_vx, partial_Omega_vy, partial_Omega_vz]
        
        # Partial derivative with respect to argument of periapsis (w)
        partial_w_x = -state_vector[1]*np.cos(i) - state_vector[2]*np.sin(i)*np.cos(Omega)
        partial_w_y = state_vector[0]*np.cos(i) - state_vector[2]*np.sin(i)*np.sin(Omega)
        partial_w_z = state_vector[0]*np.sin(i)*np.cos(Omega) + state_vector[1]*np.sin(i)*np.sin(Omega)
        partial_w_vx = -state_vector[4]*np.cos(i) - state_vector[5]*np.sin(i)*np.cos(Omega)
        partial_w_vy = state_vector[3]*np.cos(i) - state_vector[5]*np.sin(i)*np.sin(Omega)
        partial_w_vz = state_vector[3]*np.sin(i)*np.cos(Omega) + state_vector[4]*np.sin(i)*np.sin(Omega)
        
        partial_w = [partial_w_x, partial_w_y, partial_w_z, partial_w_vx, partial_w_vy, partial_w_vz]
        
        # Partial derivative with respect to mean anomaly (M)
        n = (mu/a**3)**0.5
        factor = -(mu*a**3)**0.5/r**3
        
        partial_M_x = (1/n) * state_vector[3]
        partial_M_y = (1/n) * state_vector[4]
        partial_M_z = (1/n) * state_vector[5]
        partial_M_vx = factor * state_vector[0]
        partial_M_vy = factor * state_vector[1]
        partial_M_vz = factor * state_vector[2]
        
        partial_M = [partial_M_x, partial_M_y, partial_M_z, partial_M_vx, partial_M_vy, partial_M_vz]
        
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
        Je2c[0, 0] = 1/(1-e)
        Je2c[0, 1] = q/(1-e)**2
        JX2c = np.matmul(J, Je2c)
        
        return JX2c
    
    @staticmethod
    def geo2rec(lon: float, lat: float, alt: float) -> np.ndarray:
        """
        Convert geodetic coordinates (latitude, longitude, altitude) to rectangular
        coordinates in Earth-fixed frame (ITRF93).
        
        This function converts geographic coordinates on Earth's surface to Cartesian
        coordinates in the Earth-fixed reference frame, accounting for Earth's ellipsoidal
        shape.
        
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
        numpy.ndarray
            Position vector [x, y, z] in Earth-fixed frame (ITRF93) [km].
            
        Examples
        --------
        >>> from multineas.orbit import OrbitalCoordinates
        >>> 
        >>> oc = OrbitalCoordinates()
        >>> # Chelyabinsk impact location
        >>> lon = 61.1  # degrees
        >>> lat = 54.8  # degrees
        >>> alt = 30.0  # km
        >>> 
        >>> r = oc.geo2rec(lon, lat, alt)
        >>> print(f"Position vector: {r} km")
        >>> print(f"Magnitude: {np.linalg.norm(r):.2f} km")
        
        Notes
        -----
        This function uses SPICE routines to account for Earth's ellipsoidal shape.
        The equatorial and polar radii are obtained from SPICE kernel data.
        """
        deg = np.pi / 180
        
        lon_rad = lon * deg
        lat_rad = lat * deg
        
        # Get Earth's radii from SPICE
        n, props = spy.bodvrd('399', 'RADII', 3)
        RE_spice = props[0]  # Equatorial radius
        RP_spice = props[2]  # Polar radius
        f_spice = (RE_spice - RP_spice) / RE_spice  # Flattening coefficient
        
        # Convert geodetic to rectangular coordinates
        r_earth_fixed = spy.georec(lon_rad, lat_rad, alt, RE_spice, f_spice)
        
        return r_earth_fixed
    
    @staticmethod
    def geo2eclip(lon: float, lat: float, alt: float, 
                  date: Optional[str] = None, 
                  et: Optional[float] = None,
                  frame: str = 'ITRF93') -> np.ndarray:
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
        deg = np.pi / 180
        
        lon_rad = lon * deg
        lat_rad = lat * deg
        
        # Get Earth's radii from SPICE
        n, props = spy.bodvrd('399', 'RADII', 3)
        RE_spice = props[0]
        RP_spice = props[2]
        f_spice = (RE_spice - RP_spice) / RE_spice
        
        # Convert geodetic to rectangular coordinates in Earth-fixed frame
        r_earth_fixed = spy.georec(lon_rad, lat_rad, alt, RE_spice, f_spice)
        
        # Convert ephemeris time if date is provided
        if date is not None:
            if et is not None:
                raise ValueError("Provide either 'date' or 'et', not both")
            et = spy.utc2et(date)
        elif et is None:
            raise ValueError("Either 'date' or 'et' must be provided")
        
        # Transform from Earth-fixed to Ecliptic J2000 frame
        M_ecl = spy.pxform(frame, 'ECLIPJ2000', et)
        r_earth_ecl = spy.mxv(M_ecl, r_earth_fixed)
        
        return r_earth_ecl
    
    @staticmethod
    def get_velocity_ecliptic(vx: float, vy: float, vz: float,
                              lon: float, lat: float, alt: float,
                              date: Optional[str] = None,
                              et: Optional[float] = None) -> np.ndarray:
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
        r = OrbitalCoordinates.geo2rec(lon, lat, alt)
        
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
        mx = spy.pxform('ITRF93', 'ECLIPJ2000', et)
        v_eclip = spy.mxv(mx, v_E)
        
        return v_eclip


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


def get_pre_impact_orbital_elements(lon: float, lat: float, alt: float,
                                    vx: float, vy: float, vz: float,
                                    date: str, mu: float) -> Tuple[float, float, float, 
                                                                   float, float, float, float]:
    """
    Compute pre-impact orbital elements from impact location and velocity.
    
    Given the geographic location and velocity of an impact event, this function
    computes the orbital elements of the impactor before it entered Earth's
    atmosphere. The function:
    1. Converts geographic coordinates to ecliptic J2000 coordinates
    2. Transforms the velocity to account for Earth's rotation
    3. Converts the velocity to ecliptic J2000 coordinates
    4. Computes orbital elements from the position and velocity
    
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
        UTC date and time in format 'YYYY-MM-DD HH:MM:SS'.
    mu : float
        Gravitational parameter [km³/s²] or [AU³/day²]. Note: if using km/s for
        velocity, use km³/s² for mu. If using AU/day for velocity, use AU³/day².
        
    Returns
    -------
    q : float
        Periapsis distance [same units as position].
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
        Semi-major axis [same units as position].
        
    Examples
    --------
    >>> from multineas.orbit import get_pre_impact_orbital_elements
    >>> 
    >>> # Chelyabinsk impact
    >>> lon, lat, alt = 61.1, 54.8, 30.0
    >>> vx, vy, vz = 3.5, -12.8, -6.3  # km/s
    >>> date = '2013-02-15 03:20:33'
    >>> mu = 1.32712440018e11  # km³/s² (Sun)
    >>> 
    >>> q, e, i, Omega, w, M, a = get_pre_impact_orbital_elements(
    ...     lon, lat, alt, vx, vy, vz, date, mu
    ... )
    >>> print(f"Semi-major axis: {a/149597870:.2f} AU")
    >>> print(f"Eccentricity: {e:.3f}")
    >>> print(f"Inclination: {i*180/np.pi:.1f} degrees")
    
    Notes
    -----
    This function computes the orbital elements relative to the central body (Sun).
    To get the heliocentric orbit, you need to add Earth's position and velocity
    to the impactor's position and velocity. See `OrbitIntegration` class for
    a complete integration including planetary perturbations.
    """
    oc = OrbitalCoordinates()
    
    # Convert geographic coordinates to ecliptic J2000
    r_eclip = oc.geo2eclip(lon, lat, alt, date=date)
    
    # Convert velocity to ecliptic J2000
    v_eclip = oc.get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=date)
    
    # Convert to orbital elements
    # Note: SPICE expects position in [km] and velocity in [km/s] for mu in [km³/s²]
    # Or position in [AU] and velocity in [AU/day] for mu in [AU³/day²]
    q, e, i, Omega, w, M, a = oc.transformation_x_to_e(
        r_eclip[0], r_eclip[1], r_eclip[2],
        v_eclip[0], v_eclip[1], v_eclip[2],
        mu
    )
    
    return q, e, i, Omega, w, M, a