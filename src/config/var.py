# src/config/var.py
# Created on 2025-04-06
# Author: Manish Yadav

from typing import NamedTuple

class FreeEnergyResult(NamedTuple):
    """Container for energy calculation results."""
    energy: float
    id: str

MAX_WORKERS = 11  #workers for parallel processing
BATCH_SIZE = 10

PARAM_TYPE = "Olson" # "MD" or "Mixed" or "Olson"