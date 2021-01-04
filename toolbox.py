import periodictable as pt
import numpy as np

def get_scattering_length(component, energy):
    scattering_length = 0 + 0j
    import scipy.constants as const

    cre = const.physical_constants["classical electron radius"][0]
    for key in component:
        scattering_length += np.multiply(pt.elements.symbol(key).xray.scattering_factors(energy=energy)[0], cre * component[key])
        scattering_length += (np.multiply(pt.elements.symbol(key).xray.scattering_factors(energy=energy)[1], cre * component[key]) * 1j)
    return scattering_length * 1e10
