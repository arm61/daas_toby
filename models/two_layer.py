import numpy as np
from refnx.analysis import possibly_create_parameter, Parameters
from refnx.reflect import Component


class TwoLayer(Component):
    """
    The class to describe a two-layer model for the lipid monolayer.

    Parameters
    ----------
    bs: float, array_like
        The scattering lengths for the head and tail components
    name: string
        A name for the monolayer
    """

    def __init__(self, bs, name="two_layer"):
        super(TwoLayer, self).__init__()
        if isinstance(bs[0], complex):
            self.b_real_h = possibly_create_parameter(
                bs[0].real, "{} - b_real_head".format(name)
            )
            self.b_imag_h = possibly_create_parameter(
                bs[0].imag, "{} - b_imag_head".format(name)
            )
        else:
            self.b_real_h = possibly_create_parameter(
                bs[0], "{} - b_real_head".format(name)
            )
            self.b_imag_h = possibly_create_parameter(
                0, "{} - b_imag_head".format(name)
            )
        if isinstance(bs[1], complex):
            self.b_real_t = possibly_create_parameter(
                bs[1].real, "{} - b_real_tail".format(name)
            )
            self.b_imag_t = possibly_create_parameter(
                bs[1].imag, "{} - b_imag_tail".format(name)
            )
        else:
            self.b_real_t = possibly_create_parameter(
                bs[1], "{} - b_real_tail".format(name)
            )
            self.b_imag_t = possibly_create_parameter(
                0, "{} - b_imag_tail".format(name)
            )

        self.mol_vol_h = possibly_create_parameter(
            100, "{} - molecular_volume_head".format(name)
        )
        self.mol_vol_t = possibly_create_parameter(
            100, "{} - molecular_volume_tail".format(name)
        )

        self.thick_h = possibly_create_parameter(
            100, "{} - thickness_head".format(name)
        )
        self.thick_t = possibly_create_parameter(
            100, "{} - thickness_tail".format(name)
        )

        self.phi_h = possibly_create_parameter(0.5, "{} - solvation_head".format(name))
        self.phi_t = possibly_create_parameter(0, "{} - solvation_tail".format(name))

        self.rough_h_t = possibly_create_parameter(
            3.3, "{} - roughness_head_tail".format(name)
        )
        self.rough_t_a = possibly_create_parameter(
            3.3, "{} - roughness_tail_air".format(name)
        )
        
        self.solv_sld = possibly_create_parameter(9.45, "{} - solvent sld".format(name))

        self.reverse_monolayer = False

        self.name = name

    def slabs(self, structure=None):
        """
        Returns
        -------
        slab_model = array of np.ndarray
            Slab representation of monolayer made up of two layers
        """
        layers = np.zeros((2, 5))

        layers[0, 0] = self.thick_t
        layers[0, 1] = self.b_real_t * 1e6 / self.mol_vol_t * (1-self.phi_t) + 0 * self.phi_t
        layers[0, 2] = self.b_imag_t * 1e6 / self.mol_vol_t * (1-self.phi_t) + 0 * self.phi_t
        layers[0, 3] = self.rough_t_a
        layers[0, 4] = 0

        layers[1, 0] = self.thick_h
        layers[1, 1] = self.b_real_h * 1e6 / self.mol_vol_h * (1-self.phi_h) + self.solv_sld * self.phi_h
        layers[1, 2] = self.b_imag_h * 1e6 / self.mol_vol_h * (1-self.phi_h) + 0 * self.phi_h
        layers[1, 3] = self.rough_h_t
        layers[1, 4] = 0

        if self.reverse_monolayer:
            layers = np.flipud(layers)
            layers[:, 3] = layers[::-1, 3]


        return layers

    @property
    def parameters(self):
        p = Parameters(name=self.name)
        p.extend(
            [
                self.b_real_h,
                self.b_imag_h,
                self.b_real_t,
                self.b_imag_t,
                self.thick_h,
                self.thick_t,
                self.mol_vol_h,
                self.mol_vol_t,
                self.rough_h_t,
                self.rough_t_a,
                self.phi_h,
                self.phi_t,
            ]
        )
        return p

    def logp(self):
        if self.phi_h >= 1 or self.phi_h < 0:
            return -np.inf
        if self.phi_t >= 1 or self.phi_t < 0:
            return -np.inf
        return 0

def set_constraints(lipids, structures):
    for i in range(1, len(lipids)):
        lipids[i].thick_h.constraint = lipids[0].thick_h
        lipids[i].thick_t.constraint = lipids[0].thick_t
        lipids[i].mol_vol_h.constraint = lipids[0].mol_vol_h
        lipids[i].mol_vol_t.constraint = lipids[0].mol_vol_t
        lipids[i].rough_h_t.constraint = lipids[0].rough_h_t
        lipids[i].rough_t_a.constraint = lipids[0].rough_t_a
        lipids[i].phi_h.constraint = lipids[0].phi_h
        lipids[i].phi_t.constraint = lipids[0].phi_t
    for i in range(1, len(structures)):
        structures[i][-1].rough.constraint = structures[0][-1].rough
    return lipids, structures
