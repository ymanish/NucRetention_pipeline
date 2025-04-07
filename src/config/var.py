# src/config/var.py
# Created on 2025-04-06
# Author: Manish Yadav

from typing import NamedTuple

class FreeEnergyResult(NamedTuple):
    """Container for energy calculation results."""
    energy: float
    id: str


