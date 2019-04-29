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
parser = argparse.ArgumentParser(description = "finally compute significance as fun tion of depth")
parser.add_argument('-nu_type', default="reactor_Ginna",type=str, dest = 'nu_type')
parser.add_argument('-dist',    default=10,    type=int, dest = 'dist')
args = parser.parse_args()

#detector config
num_prot = 6*9 + 1*6 #number of protons per molecule
volume=1. #1m3
rho = 999.8395 #density in kg/m3, from spec sheet
mass = volume*rho #extract total mass
molar_mass = 78.114/1000 #molar mass in kg/mol
mols = mass/molar_mass #number of mols
na = 6.0221415e23 #avogadros number in per mols
molecules = mols * na #number of molecules
total_protons = molecules * num_prot #total number of protons inside our detector

def convert_neutrinos_to_interactions_per_second(neutrinos):
    final_neutrinos = neutrinos
    for neut in final_neutrinos:
        energy = neut[0]
        sigma = 6.7e-46 * energy
        total_ints = sigma * total_protons
        neut[1] = neut[1] * total_ints
    return final_neutrinos

def convert_muons_to_interactions_per_second(muons):
    final_muons = muons
    for muon in final_muons:
        muon[1]= muon[1] * volume 
    return final_muons

def sig1(reactor,solar, muons):
    #sum over total interactions muliply by weight of each interaction
    total_muons= np.sum(muons[:,1]) * muons[:,2][0] / 4. #4 energy slices
    total_solar= np.sum(solar[:,1]) * solar[:,2][0]
    total_reactor= np.sum(reactor[:,1]) * reactor[:,2][0]
    print total_reactor, total_solar, total_muons
    return total_reactor/np.sqrt(total_solar+total_muons)
def sig2(reactor,solar, muons):
    total_muons= np.sum(muons[:,1]) * muons[:,2][0]/4. #4 energy slices
    total_solar= np.sum(solar[:,1]) * solar[:,2][0]
    total_reactor= np.sum(reactor[:,1]) * reactor[:,2][0]
    return total_reactor/np.sqrt(total_reactor+total_solar+total_muons)

t_solar=np.load("data/solar_dist_%s.npy"%args.dist)
t_reactor=np.load("data/%s_dist_%s.npy"%(args.nu_type,args.dist))
tot_solar = convert_neutrinos_to_interactions_per_second(t_solar)
tot_reactor = convert_neutrinos_to_interactions_per_second(t_reactor)
depth = [10, 100, 1000]#,10000]
sigs_1 = []
sigs_2 = []
for depths in depth:
    t_muons=np.load("data/muons_at_depth_%s.npy"%depths)
    tot_muons = convert_muons_to_interactions_per_second(t_muons)
    sig_1 = sig1(tot_reactor, tot_solar, tot_muons)
    sigs_1.append(np.log10(sig_1))
    sig_2 = sig2(tot_reactor, tot_solar, tot_muons)
    sigs_2.append(np.log10(sig_2))
    plt.plot(np.log10(depths),np.log10(sig_1),label="Depth: %s m"%depths, marker="o")
    #plt.plot(np.log10(depths),np.log10(sig_2),  marker="x")
    muons = 0

#extract optimal depths
depth_scan = np.linspace(np.log10(depth[0]),np.log10(depth[2]),1000)
params = np.polyfit(np.log10(depth),sigs_1,1)
threshold_sig = np.log10(5)
threshold_depth = (threshold_sig - params[1])/params[0]
lower_lim = 0
if(threshold_depth < np.log10(depth[0])): lower_lim = threshold_depth
else : lower_lim= np.log10(depth[0])
depth_scan = np.linspace(lower_lim,np.log10(depth[2]),1000)
sig_scan = depth_scan*params[0] +params[1]
#finish plotting
plt.plot(threshold_depth,threshold_sig,label="5sigma at %3.1fm"%10**threshold_depth,marker="*", color="red")
plt.plot(depth_scan,sig_scan,color="red")
plt.title("Signficance at Distance %s m"%args.dist)
plt.ylabel("log Significance")
plt.xlabel("log Depth of Detector [m]")
plt.legend(loc="upper left")
plt.savefig('plots/significance_at_distance_%s.pdf'%args.dist)
plt.clf()

np.save("data/optimal_depth_for_%sm.npy"%args.dist,[10**threshold_depth])
