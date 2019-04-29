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


events=np.load("data/muons_at_surface.npy")
plt.loglog(events[:,0],events[:,1],label="Muons at Surface", marker=".",linestyle='None')
plt.ylabel(r"Flux of Muons [$\frac{1}{s*m^{2}}$]")
plt.xlabel("Muon Energy [GeV]")
plt.legend(loc="upper right")
plt.savefig('plots/muon_fluxes.pdf')
plt.savefig('plots/muon_fluxes.png')
plt.clf()

events=np.load("data/muons_at_depth_10.npy")
plt.loglog(events[:,0],events[:,1],label="Muons at depth 10m", marker=".",linestyle='None')
events=np.load("data/muons_at_depth_100.npy")
plt.loglog(events[:,0],events[:,1],label="Muons at depth 100m", marker=".",linestyle='None')
events=np.load("data/muons_at_depth_1000.npy")
plt.loglog(events[:,0],events[:,1],label="Muons at depth 1km", marker=".",linestyle='None')
plt.ylabel(r"Flux of Muons [$\frac{1}{s*m^{2}}$]")
plt.xlabel("Muon Energy [GeV]")
plt.legend(loc="upper right")
plt.savefig('plots/muon_fluxes_at_depth_1km.png')
plt.savefig('plots/muon_fluxes_at_depth_1km.pdf')
plt.clf()
