"""Bremsstrahlung power and coronal impurity radiation fits.

Mavrin fits from A. A. Mavrin (2018).

Reference: https://doi.org/10.1080/10420150.2018.1462361,
Eqs. (4), (5), (7) and Tables 1--4, 6, 8, 10, 12, 14, 16, 18.
The cooling coefficients include bremsstrahlung, line and recombination
radiation. They are per electron per impurity nucleus, in W m^3.
"""

import numpy as np


# Atomic number: (internal temperature boundaries [keV], coefficient rows).
# Each row contains A0, ..., A4 for log10(Lz), in ascending powers of
# X = log10(Te [keV]). All fits span 0.1--100 keV.
MAVRIN_COEFFICIENTS = {
    2: ((), (
        (-3.5551e1, 3.1469e-1, 1.0156e-1, -9.3730e-2, 2.5020e-2),
    )),
    3: ((), (
        (-3.5115e1, 1.9475e-1, 2.5082e-1, -1.6070e-1, 3.5190e-2),
    )),
    4: ((), (
        (-3.4765e1, 3.7270e-2, 3.8363e-1, -2.1384e-1, 4.1690e-2),
    )),
    6: ((0.5,), (
        (-3.4738e1, -5.0085, -1.2788e1, -1.6637e1, -7.2904),
        (-3.4174e1, -3.6687e-1, 6.8856e-1, -2.9191e-1, 4.4470e-2),
    )),
    7: ((0.5, 2.0), (
        (-3.4065e1, -2.3614, -6.0605, -1.1570e1, -6.9621),
        (-3.3899e1, -5.9668e-1, 7.6272e-1, -1.7160e-1, 5.8770e-2),
        (-3.3913e1, -5.2628e-1, 7.0047e-1, -2.2790e-1, 2.8350e-2),
    )),
    8: ((0.3,), (
        (-3.7257e1, -1.5635e1, -1.7141e1, -5.3765, 0.0),
        (-3.3640e1, -7.6211e-1, 7.9655e-1, -2.0850e-1, 1.4360e-2),
    )),
    10: ((0.7, 5.0), (
        (-3.3132e1, 1.7309, 1.5230e1, 2.8939e1, 1.5648e1),
        (-3.3290e1, -8.7750e-1, 8.6842e-1, -3.9544e-1, 1.7244e-1),
        (-3.3410e1, -4.5345e-1, 2.9731e-1, 4.3960e-2, -2.6930e-2),
    )),
    18: ((0.6, 3.0), (
        (-3.2155e1, 6.5221, 3.0769e1, 3.9161e1, 1.5353e1),
        (-3.2530e1, 5.4490e-1, 1.5389, -7.6887, 4.9806),
        (-3.1853e1, -1.6674, 6.1339e-1, 1.7480e-1, -8.2260e-2),
    )),
    36: ((0.447, 2.364), (
        (-3.4512e1, -2.1484e1, -4.4723e1, -4.0133e1, -1.3564e1),
        (-3.1399e1, -5.0091e-1, 1.9148, -2.5865, -5.2704),
        (-2.9954e1, -6.3683, 6.6831, -2.9674, 4.8356e-1),
    )),
    54: ((0.5, 2.5, 10.0), (
        (-2.9303e1, 1.4351e1, 4.7081e1, 5.9580e1, 2.5615e1),
        (-3.1113e1, 5.9339e-1, 1.2808, -1.1628e1, 1.0748e1),
        (-2.5813e1, -2.7526e1, 4.8614e1, -3.6885e1, 1.0069e1),
        (-2.2138e1, -2.2592e1, 1.9619e1, -7.5181, 1.0858),
    )),
    74: ((1.5, 4.0), (
        (-3.0374e1, 3.8304e-1, -9.5126e-1, -1.0311, -1.0103e-1),
        (-3.0238e1, -2.9208, 2.2824e1, -6.3303e1, 5.1849e1),
        (-3.2153e1, 5.2499, -6.2740, 2.6627, -3.6759e-1),
    )),
}


def mavrin_cooling_rate(atomic_number, electron_temperature):
    """Return the total impurity cooling coefficient Lz [W m^3].

    ``electron_temperature`` is a scalar or array in eV. Temperatures
    outside 100--100000 eV use the coefficient at the nearest endpoint;
    no polynomial extrapolation is performed. At an internal boundary,
    the higher-temperature interval is used. The density-independent
    coronal fits were derived for ne = 1e20 m^-3; impurity power density
    is ne * n_imp * Lz, where n_imp includes all charge states.
    """
    if atomic_number not in MAVRIN_COEFFICIENTS:
        raise ValueError(f'No Mavrin radiation fit for atomic number {atomic_number}')
    temperature = np.asarray(electron_temperature, dtype=float)
    if not np.all(np.isfinite(temperature)):
        raise ValueError('Mavrin radiation requires finite electron temperatures')
    temperature_keV = np.clip(temperature, 100.0, 100000.0) * 1e-3
    boundaries, coefficients = MAVRIN_COEFFICIENTS[atomic_number]
    interval = np.searchsorted(boundaries, temperature_keV, side='right')
    coefficients = np.asarray(coefficients)[interval]
    x = np.log10(temperature_keV)
    log_rate = coefficients[..., 4]
    for power in range(3, -1, -1):
        log_rate = log_rate * x + coefficients[..., power]
    return 10.0**log_rate


def BremsstrahlungPower(Zi, ni, ne, Te):
    """Compute Bremsstrahlung power density for electrons slowed by an ion.

    This routine computes the Bremsstrahlung irradiated power density
    according to Freidberg's Plasma Physics and Fusion Energy,
    page 56, Eq. (3.42)

    Parameters
    ----------
    Zi : integer
        Charge number of the ion
    ni : real
        Ion density [m^-3]
    ne : real
        Electron density [m^-3]
    Te : real
        Electron temperature [eV]
    Returns
    -------
    SB : real
        Bremsstrahlung power density [W/m^3]
    """

    CB = 5.35e3
    n20 = ne / 1e20
    Tk = Te / 1e3

    Zeff = Zi*Zi * ni/ne

    SB = CB * Zeff * n20**2 * np.sqrt(Tk) # W/m^3

    return SB
