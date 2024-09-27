import pickle
# NLopt is the optimization libary dadi uses
import nlopt
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
import dadi
from dadi import Numerics, Integration, PhiManip, Spectrum

def vic_delay_size_prebottle(params, ns, pts):
    """
    nu0, Tb, s, nu1, nu2, Tpre, Tmig, m12, m21 = params
    
    Split into two vicariant populations, with asymmetrical migration after some time has passed post split. The split is preceeded by a population size change. The onset of migration coincides with a population size change.

    nu0: Size after bottleneck
    Tb: Time in the past of bottleneck before split (in units of 2*Na generations)
    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of population 1.
    nu2: Final size of population 2.
    Tpre: Time in the past after split but before migration (in units of 2*Na generations) 
    Tmig: Time in the past after migration starts (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2 (2*Na*m21)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu0, Tb, s, nu1, nu2, Tpre, Tmig, m12, m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, Tb, nu0)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    s2 = 1-s
    phi = Integration.two_pops(phi, xx, Tpre, s, s2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tmig, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

vic_delay_size_prebottle.__param_names__ = ["nu0", "Tb", "s", "nu1", "nu2", "Tpre", "Tmig", "m12", "m21"]
