# src/config/path.py
# Author: MY
# Created on: 2025-04-02


from pathlib import Path
import datetime

SRC_DIR = Path(__file__).parent.parent.parent
DATA_DIR = SRC_DIR / 'data'
BACKEND_DIR = SRC_DIR / 'backend'
RESULTS_DIR = SRC_DIR / 'output'
PARAM_DIR = SRC_DIR / 'parameters'
TEST_DIR = SRC_DIR / 'tests'
EXAMPLES_DIR = SRC_DIR / 'examples'




