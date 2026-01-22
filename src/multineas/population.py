
import numpy as np
from .util import Util

class Asteroid(object):
    """
    Class to manipulate asteroid properties.
    """
    
    # Constants for diameter-Hmag conversion
    C_std = 1329
    a_std = 0.15
    
    @staticmethod
    def diameter_to_hmag(diameter, albedo=a_std):
        """
        Convert diameter to absolute magnitude.
        
        Parameters
        ----------
        diameter : float or array
            Diameter in km.
        albedo : float, optional
            Albedo (default 0.15).
            
        Returns
        -------
        H : float or array
            Absolute magnitude.
        
        Notes
        -----
        The formula used for conversion is:
        
        $$
        D = \\frac{C}{\\sqrt{p_v}} 10^{-H/5}
        $$
        
        where $D$ is diameter in km, $p_v$ is albedo, and $H$ is absolute magnitude.
        The constant $C = 1329$ km.
        
        We can invert this to get $H$:
        
        $$
        H = -5 \\log_{10} \\left( \\frac{D \\sqrt{p_v}}{C} \\right)
        $$
        
        Examples
        --------
        >>> D = 1.0
        >>> Asteroid.diameter_to_hmag(D)
        17.7...
        
        Attr. [HC]
        """
        return -5 * np.log10(diameter * np.sqrt(albedo) / Asteroid.C_std)

    @staticmethod
    def hmag_to_diameter(H, albedo=a_std):
        """
        Convert absolute magnitude to diameter.
        
        Parameters
        ----------
        H : float or array
            Absolute magnitude.
        albedo : float, optional
            Albedo (default 0.15).
            
        Returns
        -------
        diameter : float or array
            Diameter in km.
            
        Notes
        -----
        The formula used for conversion is:
        
        $$
        D = \\frac{C}{\\sqrt{p_v}} 10^{-H/5}
        $$
        
        where $D$ is diameter in km, $p_v$ is albedo, and $H$ is absolute magnitude.
        The constant $C = 1329$ km.
        
        Examples
        --------
        >>> H = 18.0
        >>> Asteroid.hmag_to_diameter(H)
        0.89...
        
        Attr. [HC]
        """
        return (Asteroid.C_std / np.sqrt(albedo)) * 10**(-H / 5)
