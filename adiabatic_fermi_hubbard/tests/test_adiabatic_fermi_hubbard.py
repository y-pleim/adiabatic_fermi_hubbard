"""
Unit and regression test for the adiabatic_fermi_hubbard package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import adiabatic_fermi_hubbard


def test_adiabatic_fermi_hubbard_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "adiabatic_fermi_hubbard" in sys.modules
