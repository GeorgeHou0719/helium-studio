"""Original engineering conductivity anchors; optional, not NIST data.

Extracted without changing the legacy cryostat interpolation. The public
gas-property wrapper enforces its 1-300 K scope. Default EDA physics uses NIST.
"""
import numpy as np

_T = np.array([1., 2., 4., 10., 20., 50., 100., 200., 300.])
_K = np.array([.0040, .0070, .0093, .0154, .0259, .0468, .0732, .1170, .1518])


def kappa_gas(temperature):
    temperature = np.clip(np.asarray(temperature, dtype=float), _T[0], _T[-1])
    return np.exp(np.interp(np.log(temperature), np.log(_T), np.log(_K)))
