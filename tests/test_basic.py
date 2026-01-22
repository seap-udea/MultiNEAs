##################################################################
#                                                                #
# MultiNEAs: Numerical tools for near-earth asteroid dynamics   #
#            and population                                      #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)        #
##################################################################

"""
Basic tests for MultiNEAs package.

This module contains basic unit tests to verify package installation
and core functionality.
"""

import unittest
import sys
import os

# Add src directory to path for testing
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

import multineas as mn


class TestPackageImport(unittest.TestCase):
    """Test basic package import and metadata."""
    
    def test_import(self):
        """Test that the package can be imported."""
        self.assertIsNotNone(mn)
    
    def test_version(self):
        """Test that version is defined."""
        self.assertTrue(hasattr(mn, '__version__'))
        self.assertIsInstance(mn.__version__, str)
    
    def test_author(self):
        """Test that author information is defined."""
        self.assertTrue(hasattr(mn, '__author__'))
        self.assertIsInstance(mn.__author__, str)
    
    def test_email(self):
        """Test that email is defined."""
        self.assertTrue(hasattr(mn, '__email__'))
        self.assertIsInstance(mn.__email__, str)


class TestMultiNEAsBase(unittest.TestCase):
    """Test the MultiNEAsBase class."""
    
    def test_base_class_instantiation(self):
        """Test that base class can be instantiated."""
        obj = mn.MultiNEAsBase()
        self.assertIsNotNone(obj)
    
    def test_str_method(self):
        """Test string representation."""
        obj = mn.MultiNEAsBase()
        str_repr = str(obj)
        self.assertIsInstance(str_repr, str)
    
    def test_repr_method(self):
        """Test detailed representation."""
        obj = mn.MultiNEAsBase()
        repr_str = repr(obj)
        self.assertIsInstance(repr_str, str)
        self.assertIn('MultiNEAsBase', repr_str)


if __name__ == '__main__':
    unittest.main()
