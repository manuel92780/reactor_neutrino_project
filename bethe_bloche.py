'''

This problem will create Bethe-Bloch curves for different charged particles overlaid
The only energy loss processes involved here are due to ionization
Use SI units when editing
Manuel Silva - manuel.silva@wisc.edu
'''
import numpy as np
import matplotlib.pyplot as plt

#constants of nature
c  = 299792458 #[m/s]
me = 0.5109989 #[MeV/c^2]
mu = 105.65837 #[MeV/c^2] 
mp = 938.27208 #[MeV/c^2] 
ma = 3727.3793 #[MeV/c^2] 
re = 2.8179403*1e-13 #[cm]
alpha = 1/137.0

def gamma(energy,mass):
    gamm = energy/mass
    return gamm

def beta(gamm):
    bet = np.sqrt((gamm*gamm-1)/(gamm*gamm))
    return bet

#granite constans
p =      2.65 #[g/cm^3]
Z_nucl = 14+16 #[charge of Nucleus]
Ar  =     60.08 #[amu]
Ar_g =     Ar * 1.66054e-24 #[g]
ne = Z_nucl*p/Ar_g #[N/cm^3]=number density of argon atoms * number of electrons
I =      Z_nucl * 10 * 1e-6#[MeV]
plasa_energy = np.sqrt(4*np.pi*ne*re*re*re)*me/alpha#[MeV] 

def Tmax(particle, energy, gamm, bet):
    T = 0
    if("electron" in particle):
        T = energy
    else:#literally anything not an electron
        T = 2 * me * bet*bet*gamm*gamm
    return T
def delta(beta,gamma):
    delt = np.log(plasa_energy/I) + np.log(beta*gamma) - 0.5
    return delt

def dEoverdx(particle, energy, mass):
    Z = 0.0; gamm = 0.0; bet = 0.0; Tmx = 0.0;
    gamm = gamma(energy,mass)
    bet  = beta(gamm)
    Tmx  = Tmax(particle,energy,gamm,bet)
    if("alpha" in particle): Z = 2
    else:                    Z = 1
    constant = (p * Z_nucl * Z * Z / (Ar * bet*bet)) * 0.307 #[g/cm^3] * [MeV cm^2 / g] = [MeV/cm]
    first  = 0.5*np.log((2*me*bet*bet*gamm*gamm*Tmx)/(I*I))
    second = bet*bet
    third  = delta(bet,gamm)
    dedx = constant*(first-second-third)
    return dedx

def find_energies(energies, stopping_power):
    min_dedx = np.min(stopping_power)
    range_dedx = 1.5 * min_dedx
    energy_point = []
    stopping_point = []
    for i in range(1,len(energies)):
        if(stopping_power[i] < range_dedx) and (stopping_power[i-1] > range_dedx): 
            energy_point.append(energies[i]); stopping_point.append(stopping_power[i]);
        if(stopping_power[i] > range_dedx) and (stopping_power[i-1] < range_dedx):
            energy_point.append(energies[i]); stopping_point.append(stopping_power[i]);
    if len(energy_point) < 2: 
        energy_point.append(energy_point[0])
        stopping_point.append(stopping_point[0])
        #energy_point = np.round(energy_point)
    return energy_point, stopping_point

#get the loops setup
particles = ["electrons", "muons", "protons", "alphas"]
masses    = [me, mu, mp, ma] #MeV/c^2 for all particles
N_p = 0 #particle number
emin = 0.1   #[MeV]
emax = 1e8 #[MeV]
ebins = 100000
logebins = np.logspace(np.log10(emin),np.log10(emax), ebins)
for N_p in range(0,len(particles)-2):
    real_e = []
    de_dx = []
    for e in logebins:
        if(e < masses[N_p]): continue
        dedx = dEoverdx(particles[N_p], e, masses[N_p])
        real_e.append(e)
        de_dx.append(dedx)
    #plot the result of a single particle
    energies, stopping = find_energies(real_e, de_dx)
    if("electron" in particles[N_p]) or ("muon" in particles[N_p]):
        label = particles[N_p]#+", E=["+str(round(energies[0],2))+"MeV,"+str(round(energies[1]/1000.0,1))+"GeV]"
    else:
        if(energies[0] > 1000):
            label = particles[N_p]+", E=["+str(round(energies[0]/1000.0,1))+","+str(round(energies[1]/1000.0,1))+"]GeV"
        else:
            label = particles[N_p]+", E=["+str(round(energies[0],1))+","+str(round(energies[1],1))+"]MeV"
        
    plt.plot(real_e, de_dx,label=label)
    #plt.scatter(energies, stopping)

plt.title('Stopping power in Granite')
plt.ylabel('dE/dx [MeV/cm]')
plt.xlabel('Total Energy of Particle [MeV]')
plt.legend(loc='upper right')
plt.xscale('log')
plt.ylim(0, 35)
plt.grid(linestyle=':', which='both', linewidth=0.5)
plt.savefig("Bethe-Bloch_various_particles.pdf")
plt.clf()
