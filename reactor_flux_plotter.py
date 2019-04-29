import numpy as np
import glob
import matplotlib
matplotlib.use('agg')
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import colors
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter


import argparse
parser = argparse.ArgumentParser(description = "inject neutrinos into our target medium at certain distances")
parser.add_argument('-nu_type', default="reactor_Ginna",type=str, dest = 'nu_type')
parser.add_argument('-depth',   default=1e3,    type=int, dest = 'depth')
args = parser.parse_args()

for dists in [10,100,1000]:
    events=np.load("data/%s_dist_%s.npy"%(args.nu_type,dists ))
    plt.loglog(events[:,0],events[:,1],label="%s Neutrinos at distance %0.1fkm"%(args.nu_type,dists/1000.), marker=".",linestyle='None')
    plt.ylabel(r"Flux at depth and distance [$\frac{1}{s*m^{2}}$]")
    plt.xlabel("Neutrino Energy [MeV]")

plt.legend(loc="lower left")
plt.savefig('plots/%s_fluxes.pdf'%args.nu_type)
plt.savefig('plots/%s_fluxes.png'%args.nu_type)
