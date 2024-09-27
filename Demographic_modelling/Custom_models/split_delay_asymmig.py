import pickle
# NLopt is the optimization libary dadi uses
import nlopt
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
import dadi
from dadi import Numerics, Integration, PhiManip, Spectrum

def split_delay_asymmig(params, ns, pts):
    """
    params = (nu1,nu2,Tpre,Tmig,m12,m21)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration after some time has passed post split.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Tpre: Time in the past after split but before migration (in units of 2*Na generations) 
    Tmig: Time in the past after migration starts (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2 (2*Na*m21)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,Tpre,Tmig,m12,m21 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Tpre, nu1, nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tmig, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

split_delay_asymmig.__param_names__ = ["nu1", "nu2", "Tpre", "Tmig", "m12", "m21"]