##################################################################
#                                                                #
# MultiNEAs: Numerical tools for near-earth asteroid dynamics   #
#            and population                                      #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)        #
# Authors: Jorge I. Zuluaga, Juanita A. Agudelo                 #
# Contact: jorge.zuluaga@udea.edu.co                            #
##################################################################

"""
MultiNEAs: Numerical Tools for Near-Earth Asteroid Dynamics and Population

This package provides a comprehensive suite of tools for studying the dynamics 
and population of Near-Earth Asteroids (NEAs).

Main Features:
--------------
- Orbital dynamics calculations
- Population statistics and analysis
- Visualization tools
- Data management utilities

Usage:
------
    >>> import multineas as mn
    >>> # Your NEA analysis code here

For more information, visit: https://github.com/seap-udea/MultiNEAs
"""

import warnings
warnings.filterwarnings('ignore')

# Package metadata
__version__ = '0.3.4'
__author__ = 'Jorge I. Zuluaga, Juanita A. Agudelo'
__email__ = 'jorge.zuluaga@udea.edu.co, juanita.agudelo@udea.edu.co'
__license__ = 'AGPL-3.0-only'
__url__ = 'https://github.com/seap-udea/MultiNEAs'

# Version information
VERSION = __version__
AUTHOR = __author__
EMAIL = __email__

import os
ROOTDIR = os.path.dirname(os.path.abspath(__file__))

class MultiNEAsBase:
    """
    Base class for MultiNEAs package.
    
    All major classes in the package inherit from this base class,
    providing common functionality and attributes.

    Attribution
    -----------
    [VC] This element was vibe coded using ChatGPT 5.2 in Antigravity.
    """
    def __init__(self):
        pass
    
    def __str__(self):
        """String representation of the object."""
        return str({k: v for k, v in self.__dict__.items() if not k.startswith('_')})
    
    def __repr__(self):
        """Detailed representation of the object."""
        return f"{self.__class__.__name__}({self.__dict__})"


# Placeholder for future imports
# When you add modules, uncomment and modify these lines:
# from multineas.core import *
# from multineas.dynamics import *
# from multineas.population import *
# from multineas.visualization import *

# Package initialization message (optional, can be removed in production)
def _welcome_message():
    """Display welcome message on import (for development)."""
    print(f"Welcome to MultiNEAs v{__version__}")

_welcome_message()

# Clean up namespace
del warnings
del _welcome_message
