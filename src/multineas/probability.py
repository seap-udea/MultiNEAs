"""
Probability density functions for orbital elements and phase space.

This module provides classes and functions for computing probability density functions
(PDFs) in orbital element space and phase space (Cartesian coordinates), as well as
integrating these PDFs over various regions.
"""

import numpy as np
from typing import Union, Tuple, List
from .orbit import OrbitalCoordinates
from .multimin import FitCMND
from .util import Util


class OrbitElementsPDF:
    """
    Probability density function for orbital elements (q, e, i, Omega, w, M).
    
    This class computes the probability density P(q, e, i, Omega, w, M) for
    Near-Earth Asteroids (NEAs) based on a Composed Multivariate Normal Distribution
    (CMND) fit. The PDF is computed in the orbital element space, where the distribution
    over (q, e, i) is modeled by the CMND, and (Omega, w, M) are assumed to be
    uniformly distributed.
    
    Attributes
    ----------
    fit_cmnd : FitCMND
        Fitted CMND model for orbital elements (q, e, i).
    scales : list or numpy.ndarray
        Scales for transforming orbital elements to unbound space [q_scale, e_scale, i_scale].
        Default is [1.30, 1.00, np.pi].
    
    Examples
    --------
    >>> from multineas.probability import OrbitElementsPDF
    >>> from multineas.multimin import FitCMND
    >>> import numpy as np
    >>> 
    >>> # Load a fitted CMND model
    >>> F = FitCMND("path/to/fit.pkl")
    >>> 
    >>> # Create PDF instance
    >>> pdf = OrbitElementsPDF(F)
    >>> 
    >>> # Compute PDF for specific orbital elements
    >>> q, e, i = 0.5, 0.3, np.pi/4
    >>> prob = pdf.compute_p_e_cmnd(q, e, i)
    >>> print(f"PDF value: {prob:.6e}")
    
    Notes
    -----
    The PDF is computed as:
    P_E_CMND = P_qei(q, e, i) * P_Omega_w_M
    
    where:
    - P_qei is the CMND PDF evaluated at the transformed (q, e, i)
    - P_Omega_w_M = 1/(2π)^3 is the uniform distribution over angles
    """
    
    def __init__(self, fit_cmnd: FitCMND, scales: Union[List[float], np.ndarray] = None):
        """
        Initialize OrbitElementsPDF.
        
        Parameters
        ----------
        fit_cmnd : FitCMND
            Fitted CMND model for orbital elements.
        scales : list or numpy.ndarray, optional
            Scales for transformation [q_scale, e_scale, i_scale].
            Default is [1.35, 1.00, np.pi].
        """
        self.fit_cmnd = fit_cmnd
        if scales is None:
            self.scales = [1.30, 1.00, np.pi]
        else:
            self.scales = np.asarray(scales)
    
    def compute_p_e_cmnd(self, q: float, e: float, i: float) -> float:
        """
        Compute probability density P_E_CMND(q, e, i, Omega, w, M).
        
        The PDF is computed as the product of:
        1. P_qei(q, e, i): CMND PDF for orbital elements (q, e, i)
        2. P_Omega_w_M: Uniform distribution over angles (Omega, w, M)
        
        Parameters
        ----------
        q : float
            Periapsis distance (AU).
        e : float
            Eccentricity (dimensionless).
        i : float
            Inclination (radians).
        
        Returns
        -------
        float
            Probability density value P_E_CMND.
        
        Examples
        --------
        >>> from multineas.probability import OrbitElementsPDF
        >>> from multineas.multimin import FitCMND
        >>> import numpy as np
        >>> 
        >>> # Load fitted model
        >>> F = FitCMND("path/to/fit.pkl")
        >>> pdf = OrbitElementsPDF(F)
        >>> 
        >>> # Compute PDF for a circular orbit at 1 AU
        >>> q = 1.0
        >>> e = 0.0
        >>> i = np.pi / 6  # 30 degrees
        >>> prob = pdf.compute_p_e_cmnd(q, e, i)
        >>> print(f"PDF: {prob:.6e}")

        Attr. [HC]
        """
        # Transform orbital elements to unbound space
        element = np.array([q, e, i])
        u_element = Util.t_if(element, self.scales, Util.f2u)
        
        # Compute CMND PDF for (q, e, i)
        P_qei = self.fit_cmnd.cmnd.pdf(u_element)
        
        # Uniform distribution over angles (Omega, w, M)
        # Range: [0, 2π] for each angle
        max_angle = 2 * np.pi
        min_angle = 0
        P_WwM = 1 / (max_angle - min_angle)**3
        
        return P_qei * P_WwM


class PhaseSpacePDF:
    """
    Probability density function for phase space (Cartesian coordinates).
    
    This class computes the probability density P_X_CMND(x, y, z, vx, vy, vz) in
    Cartesian phase space by transforming to orbital elements and applying the
    appropriate Jacobian transformations.
    
    Attributes
    ----------
    orbit_coords : OrbitalCoordinates
        Instance for coordinate transformations.
    orbit_elements_pdf : OrbitElementsPDF
        PDF for orbital elements.
    q_max : float
        Maximum perihelion distance for NEA classification (AU).
    e_max : float
        Maximum eccentricity for valid orbits.
    i_max : float
        Maximum inclination (radians).
    mu : float
        Gravitational parameter (AU^3/day^2).
    
    Examples
    --------
    >>> from multineas.probability import PhaseSpacePDF, OrbitElementsPDF
    >>> from multineas.orbit import OrbitalCoordinates
    >>> from multineas.multimin import FitCMND
    >>> import numpy as np
    >>> 
    >>> # Setup
    >>> F = FitCMND("path/to/fit.pkl")
    >>> oc = OrbitalCoordinates()
    >>> oepdf = OrbitElementsPDF(F)
    >>> 
    >>> # Create phase space PDF
    >>> mu = 0.01720209895**2
    >>> pdf = PhaseSpacePDF(oc, oepdf, q_max=1.35, e_max=1.0, i_max=np.pi, mu=mu)
    >>> 
    >>> # Compute PDF at a point in phase space
    >>> x, y, z = 1.0, 0.0, 0.0
    >>> vx, vy, vz = 0.0, 6.28, 0.0
    >>> prob = pdf.compute_p_x_cmnd(x, y, z, vx, vy, vz)
    >>> print(f"PDF: {prob:.6e}")

    Attr. [HC]
    """
    
    # Constants for NEA orbit validation
    Q_MAX_NEA = 1.3  # Maximum perihelion distance for NEA classification (AU)
    E_MAX_VALID = 1.0  # Maximum eccentricity for valid elliptical orbits
    
    def __init__(self, orbit_coords: OrbitalCoordinates, orbit_elements_pdf: OrbitElementsPDF,
                 q_max: float, e_max: float, i_max: float, mu: float):
        """
        Initialize PhaseSpacePDF.
        
        Parameters
        ----------
        orbit_coords : OrbitalCoordinates
            Instance for coordinate transformations.
        orbit_elements_pdf : OrbitElementsPDF
            PDF for orbital elements.
        q_max : float
            Maximum perihelion distance (AU).
        e_max : float
            Maximum eccentricity.
        i_max : float
            Maximum inclination (radians).
        mu : float
            Gravitational parameter (AU^3/day^2).
        """
        self.orbit_coords = orbit_coords
        self.orbit_elements_pdf = orbit_elements_pdf
        self.q_max = q_max
        self.e_max = e_max
        self.i_max = i_max
        self.mu = mu
    
    def compute_p_x_cmnd(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray],
                        z: Union[float, np.ndarray], vx: Union[float, np.ndarray],
                        vy: Union[float, np.ndarray], vz: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute probability density P_X_CMND(x, y, z, vx, vy, vz).
        
        This function transforms Cartesian coordinates to orbital elements, computes
        the PDF in orbital element space, and applies the Jacobian transformation
        to get the PDF in phase space.
        
        Parameters
        ----------
        x, y, z : float or numpy.ndarray
            Position coordinates (AU).
        vx, vy, vz : float or numpy.ndarray
            Velocity coordinates (AU/day).
        
        Returns
        -------
        float or numpy.ndarray
            Probability density values, same shape as broadcasted input arrays.
        
        Examples
        --------
        >>> from multineas.probability import PhaseSpacePDF, OrbitElementsPDF
        >>> from multineas.orbit import OrbitalCoordinates
        >>> from multineas.multimin import FitCMND
        >>> import numpy as np
        >>> 
        >>> # Setup
        >>> F = FitCMND("path/to/fit.pkl")
        >>> oc = OrbitalCoordinates()
        >>> oepdf = OrbitElementsPDF(F)
        >>> mu = 0.01720209895**2
        >>> pdf = PhaseSpacePDF(oc, oepdf, q_max=1.35, e_max=1.0, i_max=np.pi, mu=mu)
        >>> 
        >>> # Single point
        >>> prob = pdf.compute_p_x_cmnd(1.0, 0.0, 0.0, 0.0, 6.28, 0.0)
        >>> 
        >>> # Array of points
        >>> x_vals = np.linspace(0.5, 1.5, 10)
        >>> y_vals = np.zeros(10)
        >>> z_vals = np.zeros(10)
        >>> vx_vals = np.zeros(10)
        >>> vy_vals = np.full(10, 6.28)
        >>> vz_vals = np.zeros(10)
        >>> probs = pdf.compute_p_x_cmnd(x_vals, y_vals, z_vals, vx_vals, vy_vals, vz_vals)

        Attr. [HC]
        """
        # Convert inputs to arrays and get broadcast shape
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        vx = np.asarray(vx)
        vy = np.asarray(vy)
        vz = np.asarray(vz)
        
        shape = np.broadcast(x, y, z, vx, vy, vz).shape
        P = np.empty(shape, dtype=float)
        
        # Flatten arrays for iteration
        x_flat = x.ravel()
        y_flat = y.ravel()
        z_flat = z.ravel()
        vx_flat = vx.ravel()
        vy_flat = vy.ravel()
        vz_flat = vz.ravel()
        
        # Process each point in the flattened arrays
        for idx in range(x_flat.size):
            # Transform from Cartesian to orbital elements
            q, e, i, Omega, w, M, a = self.orbit_coords.transformation_x_to_e(
                x_flat[idx], y_flat[idx], z_flat[idx],
                vx_flat[idx], vy_flat[idx], vz_flat[idx], self.mu
            )
            
            # Validate orbit: must be NEA (q <= Q_MAX_NEA) and elliptical (e < E_MAX_VALID)
            if q > self.Q_MAX_NEA or e >= self.E_MAX_VALID:
                P.flat[idx] = 0.0
                continue
            
            # Compute Jacobian matrices
            jacobian_x_to_e = self.orbit_coords.compute_jacobian_x_to_e(
                a, e, i, Omega, w, M, self.mu
            )
            jacobian_qei_to_QEI = compute_jacobian_qei_to_QEI(
                q, e, i, self.q_max, self.e_max, self.i_max
            )
            
            # Compute determinants
            det_X_to_E = np.linalg.det(jacobian_x_to_e)
            det_QEI = (jacobian_qei_to_QEI[0, 0] * 
                      jacobian_qei_to_QEI[1, 1] * 
                      jacobian_qei_to_QEI[2, 2])
            
            # Check for computational inconsistencies
            if det_X_to_E == 0 or not np.isfinite(det_X_to_E):
                P.flat[idx] = np.nan
                continue
            
            # Compute probability density
            inv_det_X_to_E = 1.0 / det_X_to_E
            pdf_orbital_elements = self.orbit_elements_pdf.compute_p_e_cmnd(q, e, i)
            P.flat[idx] = pdf_orbital_elements * abs(inv_det_X_to_E) * abs(det_QEI)
        
        return P.reshape(shape)


def compute_jacobian_qei_to_QEI(q: float, e: float, i: float,
                                q_max: float, e_max: float, i_max: float) -> np.ndarray:
    """
    Compute Jacobian matrix for transformation from (q, e, i) to (Q, E, I).
    
    This Jacobian corresponds to the transformation from finite intervals [0, max]
    to unbound space using the u2f transformation from Util. The transformation
    is: Q = q_max/(1 + exp(-u_q)), where u_q is the unbound variable.
    
    Parameters
    ----------
    q : float
        Periapsis distance (AU).
    e : float
        Eccentricity (dimensionless).
    i : float
        Inclination (radians).
    q_max : float
        Maximum perihelion distance (AU).
    e_max : float
        Maximum eccentricity.
    i_max : float
        Maximum inclination (radians).
    
    Returns
    -------
    numpy.ndarray
        Jacobian matrix of shape (3, 3).
    
    Examples
    --------
    >>> from multineas.probability import compute_jacobian_qei_to_QEI
    >>> import numpy as np
    >>> 
    >>> q, e, i = 0.5, 0.3, np.pi/4
    >>> q_max, e_max, i_max = 1.35, 1.0, np.pi
    >>> J = compute_jacobian_qei_to_QEI(q, e, i, q_max, e_max, i_max)
    >>> print(f"Jacobian determinant: {np.linalg.det(J):.6e}")
    
    Attr. [HC]
    """
    partialQ_q = q_max / (q * (q_max - q))
    partialQ_e = 0
    partialQ_i = 0
    partialE_q = 0
    partialE_e = e_max / (e * (e_max - e))
    partialE_i = 0
    partialI_q = 0
    partialI_e = 0
    partialI_i = i_max / (i * (i_max - i))
    
    Jacobian = np.array([
        [partialQ_q, partialQ_e, partialQ_i],
        [partialE_q, partialE_e, partialE_i],
        [partialI_q, partialI_e, partialI_i]
    ])
    return Jacobian


def integrate(center: Union[Tuple, List, np.ndarray], widths: Union[Tuple, List, np.ndarray],
              max_elements: Union[Tuple, List], pdf: PhaseSpacePDF,
              n_points: int = 8) -> float:
    """
    Compute hypercube surface integral of P_X_CMND over a 6D hypercube.
    
    This function integrates the phase space PDF over a hypercube centered at
    (x, y, z, vx, vy, vz) with specified widths using Gauss-Legendre quadrature.
    
    Parameters
    ----------
    center : tuple, list, or numpy.ndarray
        Center of hypercube (x, y, z, vx, vy, vz).
    widths : tuple, list, or numpy.ndarray
        Side lengths of hypercube (dx, dy, dz, dvx, dvy, dvz).
    max_elements : tuple, list, or numpy.ndarray
        Maximum orbital elements (q_max, e_max, i_max).
    pdf : PhaseSpacePDF
        Phase space PDF instance.
    n_points : int, optional
        Number of quadrature points per dimension (default 8).
    
    Returns
    -------
    float
        Integral value.
    
    Examples
    --------
    >>> from multineas.probability import integrate, PhaseSpacePDF, OrbitElementsPDF
    >>> from multineas.orbit import OrbitalCoordinates
    >>> from multineas.multimin import FitCMND
    >>> import numpy as np
    >>> 
    >>> # Setup
    >>> F = FitCMND("path/to/fit.pkl")
    >>> oc = OrbitalCoordinates()
    >>> oepdf = OrbitElementsPDF(F)
    >>> mu = 0.01720209895**2
    >>> pdf = PhaseSpacePDF(oc, oepdf, q_max=1.35, e_max=1.0, i_max=np.pi, mu=mu)
    >>> 
    >>> # Integrate over hypercube
    >>> center = [1.0, 0.0, 0.0, 0.0, 6.28, 0.0]
    >>> widths = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    >>> max_elements = [1.35, 1.0, np.pi]
    >>> integral = integrate(center, widths, max_elements, pdf, n_points=6)
    >>> print(f"Integral: {integral:.6e}")

    Attr. [MC]
    """
    from numpy.polynomial.legendre import leggauss
    
    x0, y0, z0, vx0, vy0, vz0 = center
    dx, dy, dz, dvx, dvy, dvz = widths
    
    # Get Gauss-Legendre points and weights for [-1, 1]
    pts, wts = leggauss(n_points)
    
    # Map points from [-1, 1] to [center-width/2, center+width/2] for each dimension
    x_pts = x0 + 0.5 * dx * pts
    y_pts = y0 + 0.5 * dy * pts
    z_pts = z0 + 0.5 * dz * pts
    vx_pts = vx0 + 0.5 * dvx * pts
    vy_pts = vy0 + 0.5 * dvy * pts
    vz_pts = vz0 + 0.5 * dvz * pts
    
    # Create meshgrid of all quadrature points
    X, Y, Z, VX, VY, VZ = np.meshgrid(
        x_pts, y_pts, z_pts, vx_pts, vy_pts, vz_pts, indexing='ij'
    )
    WX, WY, WZ, WVX, WVY, WVZ = np.meshgrid(
        wts, wts, wts, wts, wts, wts, indexing='ij'
    )
    
    # Flatten for vectorized evaluation
    Xf = X.ravel()
    Yf = Y.ravel()
    Zf = Z.ravel()
    VXf = VX.ravel()
    VYf = VY.ravel()
    VZf = VZ.ravel()
    WF = (WX * WY * WZ * WVX * WVY * WVZ).ravel()
    
    # Evaluate P at all points
    Pf = pdf.compute_p_x_cmnd(Xf, Yf, Zf, VXf, VYf, VZf)
    
    # Integral is sum(P * weight) * volume factor
    integral = np.sum(Pf * WF) * (0.5 * dx) * (0.5 * dy) * (0.5 * dz) * (0.5 * dvx) * (0.5 * dvy) * (0.5 * dvz)
    return integral


def marginalize(center: Union[Tuple, List, np.ndarray], widths: Union[Tuple, List, np.ndarray],
                fixed_values: Union[Tuple, List, np.ndarray], max_elements: Union[Tuple, List],
                pdf: PhaseSpacePDF, integrate_over: str = 'velocity', n_points: int = 8) -> float:
    """
    Compute box surface integral of P_X_CMND over a 3D box.
    
    This function integrates the phase space PDF over a 3D box in either position
    or velocity space, while keeping the other fixed. Uses Gauss-Legendre quadrature.
    
    Parameters
    ----------
    center : tuple, list, or numpy.ndarray
        Center of integration box (3D: c1, c2, c3).
    widths : tuple, list, or numpy.ndarray
        Side lengths of integration box (3D: d1, d2, d3).
    fixed_values : tuple, list, or numpy.ndarray
        Fixed values for non-integrated dimensions (3D: f1, f2, f3).
    max_elements : tuple, list, or numpy.ndarray
        Maximum orbital elements (q_max, e_max, i_max).
    pdf : PhaseSpacePDF
        Phase space PDF instance.
    integrate_over : str, optional
        'position' or 'velocity' - specifies which dimensions to integrate over (default 'velocity').
    n_points : int, optional
        Number of quadrature points per dimension (default 8).
    
    Returns
    -------
    float
        Integral value.
    
    Examples
    --------
    >>> from multineas.probability import marginalize, PhaseSpacePDF, OrbitElementsPDF
    >>> from multineas.orbit import OrbitalCoordinates
    >>> from multineas.multimin import FitCMND
    >>> import numpy as np
    >>> 
    >>> # Setup
    >>> F = FitCMND("path/to/fit.pkl")
    >>> oc = OrbitalCoordinates()
    >>> oepdf = OrbitElementsPDF(F)
    >>> mu = 0.01720209895**2
    >>> pdf = PhaseSpacePDF(oc, oepdf, q_max=1.35, e_max=1.0, i_max=np.pi, mu=mu)
    >>> 
    >>> # Integrate over velocity box with fixed position
    >>> center_v = [0.0, 6.28, 0.0]
    >>> widths_v = [0.5, 0.5, 0.5]
    >>> fixed_pos = [1.0, 0.0, 0.0]
    >>> max_elements = [1.35, 1.0, np.pi]
    >>> integral = marginalize(center_v, widths_v, fixed_pos, max_elements, pdf, 
    ...                        integrate_over='velocity', n_points=6)
    >>> print(f"Integral: {integral:.6e}")
    
    Attr. [MC]
    """
    from numpy.polynomial.legendre import leggauss
    
    if integrate_over not in ['position', 'velocity']:
        raise ValueError("integrate_over must be 'position' or 'velocity'")
    
    c1, c2, c3 = center
    d1, d2, d3 = widths
    f1, f2, f3 = fixed_values
    
    # Get Gauss-Legendre points and weights for [-1, 1]
    pts, wts = leggauss(n_points)
    
    # Map points from [-1, 1] to [center-width/2, center+width/2] for integration dimensions
    dim1_pts = c1 + 0.5 * d1 * pts
    dim2_pts = c2 + 0.5 * d2 * pts
    dim3_pts = c3 + 0.5 * d3 * pts
    
    # Create meshgrid for integration dimensions
    DIM1, DIM2, DIM3 = np.meshgrid(dim1_pts, dim2_pts, dim3_pts, indexing='ij')
    W1, W2, W3 = np.meshgrid(wts, wts, wts, indexing='ij')
    
    # Flatten for vectorized evaluation
    DIM1f = DIM1.ravel()
    DIM2f = DIM2.ravel()
    DIM3f = DIM3.ravel()
    WF = (W1 * W2 * W3).ravel()
    
    # Create constant arrays for fixed dimensions (same size as integration arrays)
    F1f = np.full_like(DIM1f, f1)
    F2f = np.full_like(DIM2f, f2)
    F3f = np.full_like(DIM3f, f3)
    
    # Assign to position and velocity based on integrate_over parameter
    if integrate_over == 'velocity':
        # Integrating over velocities, positions are fixed
        Xf, Yf, Zf = F1f, F2f, F3f  # fixed positions
        VXf, VYf, VZf = DIM1f, DIM2f, DIM3f  # integrated velocities
    else:  # integrate_over == 'position'
        # Integrating over positions, velocities are fixed
        Xf, Yf, Zf = DIM1f, DIM2f, DIM3f  # integrated positions
        VXf, VYf, VZf = F1f, F2f, F3f  # fixed velocities
    
    # Evaluate P at all points
    Pf = pdf.compute_p_x_cmnd(Xf, Yf, Zf, VXf, VYf, VZf)
    
    # Integral is sum(P * weight) * volume factor
    integral = np.sum(Pf * WF) * (0.5 * d1) * (0.5 * d2) * (0.5 * d3)
    return integral
