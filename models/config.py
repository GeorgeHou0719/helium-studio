"""Portable local runtime paths, configurable before starting the application."""
import os
from pathlib import Path
import shutil

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get('HELIUM_STUDIO_DATA_DIR', PROJECT_ROOT / 'var')).expanduser().resolve()
COMSOL_TEMPLATE = Path(os.environ.get(
    'HELIUM_STUDIO_COMSOL_TEMPLATE', DATA_DIR / 'comsol' / 'helium_ss304_pipe.mph'
)).expanduser().resolve()
_default_executable = shutil.which('comsol') or (
    'C:/Program Files/COMSOL/COMSOL61/Multiphysics/bin/win64/comsol.exe'
    if os.name == 'nt' else '/usr/local/comsol61/multiphysics/bin/comsol'
)
COMSOL_EXECUTABLE = Path(os.environ.get('COMSOL_EXECUTABLE', _default_executable)).expanduser()
