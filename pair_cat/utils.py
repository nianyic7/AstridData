c_mks = 3e8
msun_mks = 2e30
s_to_year = 3.17098e-8
year_to_s = 1.0 / s_to_year
lsun_ergs = 3.9e33
mdot_msun_yr = 1e10 / 980 / 1e6


def calc_lx(mdot):
    """
    input: mdot in Msun/yr
    output: Lx in ergs
    """
    lbol = 0.1 * mdot * msun_mks / year_to_s * c_mks**2
    lbol_lsun = lbol / 3.9e26
    k = 10.83 * (lbol_lsun / 1e10) ** 0.28 + 6.08 * (lbol_lsun / 1e10) ** (-0.02)
    return lbol / k * 1e7


def calc_lbol(mdot):
    """
    input: mdot in Msun/yr
    output: Lx in ergs
    """
    lbol = 0.1 * mdot * msun_mks / year_to_s * c_mks**2
    lbol_ergs = lbol * 1e7
    return lbol_ergs


def edd_ratio(mass, lum):
    return lum / (1.26e38 * mass)