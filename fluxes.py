import numpy as np
import bethe_bloche as bb

#neutrino constants
nu_energy_high=2.5 #10 MeV
nu_energy_low =0.5  #1 MeV
reactors = {"Ginna":580,"Grand_Guld":1500}#power output in MegaWatts eg MJ/s
MJ_to_MeV = 6.242e18 #convert 1 MW to 1 MeV/s
def gaussian(energy, mu=1.5, sig=0.5):
    return np.exp(-np.power(energy - mu, 2.) / (2 * np.power(sig, 2.)))

#muon constants
mu = 105.65837 #[MeV/c^2]
e1 = 1.000;
e2 = 927.7;
e3 = 1588.;
e4 = 416.2*1000;
e5 = 1000.*1000;
def flux_params(energy):
    if energy < e2:
        return {"C":2.950e-3,"p1":0.3061,"p2":1.2743,"p3":-0.2630,"p4":0.0252}
    if energy < e3:
        return {"C":1.781e-2,"p1":1.7910,"p2":0.3040,"p3":0.0000,"p4":0.0000}
    if energy < e4:
        return {"C":1.435e1, "p1":3.6720,"p2":0.0000,"p3":0.0000,"p4":0.0000}
    if energy < e5:
        return {"C":1.000e3, "p1":4     ,"p2":0.0000,"p3":0.0000,"p4":0.0000}
def muon_flux(energy):
    ps = flux_params(energy)
    C = ps["C"];
    p0 = ps["p1"]; p1 = ps["p2"]; p2 = ps["p3"]; p3 = ps["p4"];
    phi_temp = C * energy**(-(p0+p1*np.log10(energy)+p2*np.log10(energy)*np.log10(energy)+p3*np.log10(energy)*np.log10(energy)*np.log10(energy)))
    phi = phi_temp * 4*np.pi*energy #convert to per cm^2 per second multiply by energy and 4pi steradians
    phi = phi * 100*100 #convert to per m^2 per second
    return phi

#solar neutrino fluxes
def solar_flux(energy):
    phi = 5.05 * 10**6 #per cm^2 per second 
    phi = phi * 100*100 #convert to per m^2 per second
    phi = phi /6. #only care about electron anti-neutrinos
    phi = phi/ 5 #only using 1/5 of the energy ranges
    return phi

#reactor neutrinos fluxes
def reactor_rate(nu_type,energy):
    reactor_power = reactors[nu_type.replace("reactor_","")]
    reactor_power *= MJ_to_MeV
    n_nu_per_fission = 6. * gaussian(energy=energy)#average number neutrinos per second
    energy_per_fission = 187 #average energy per fission
    dN_dt =  reactor_power * n_nu_per_fission / energy_per_fission #neutrino rates as function of energy
    return dN_dt

def reactor_flux(nu_type,energy, distance):
    phi = reactor_rate(nu_type,energy)/(4*np.pi*distance*distance) #convert to per m^2 per second
    return phi

def simulate_neutrinos(nu_type,n_events, distance):
    events = []
    if"solar" in nu_type:
        energies = np.random.uniform(nu_energy_low,nu_energy_high,n_events )
    if "reactor" in nu_type:
        energies = np.random.uniform(0.5, 4, n_events)
    for event_i in range(len(energies)):
        energy = energies[event_i]
        weight = 1./len(energies)
        dist   = distance
        if "solar" in nu_type:
            flux   = solar_flux(energy)
        if "reactor" in nu_type:
            flux   = reactor_flux(nu_type,energy,distance)
        event =  [energy,flux,weight,dist] 
        events.append(event)
    return events

def simulate_muons(n_events, depth):
    #inject only in GeV, energy ranges
    energy_1 = np.random.uniform(e1,e2,n_events)
    energy_2 = np.random.uniform(e2,e3,n_events)
    energy_3 = np.random.uniform(e3,e4,n_events)
    energy_4 = np.random.uniform(e4,e5,n_events)
    energies = np.concatenate((energy_1,energy_2,energy_3,energy_4))
    events = []
    for event_i in range(len(energies)):
        energy = energies[event_i]
        flux   = muon_flux(energy)
        weight = 1./len(energies)
        dep    = depth*0
        flux   = muon_flux(energy)
        weight = 1./n_events
        event =  [energy,flux,weight,dep]
        events.append(event)
    return events

def propagate_muons(muons_surface,depth):
    final_muons = muons_surface
    step_dist = 10 #take steps of 10m
    steps = depth/step_dist #number of steps to run over energy losses
    for step in range(int(steps)):
        for muon in final_muons:
            energy = muon[0]
            if(energy < (1.2/1000.)): #if emergy less than 1.2MeV, end the event 
                muon[0] = 0; muon[2] = 0; 
                continue #set energy and flux to zero manually
            losses = bb.dEoverdx("muons",(energy*1000.0)+mu,mu) #energy loss in MeV/cm
            energy_lost = (losses * step_dist*100)/1000.0 #energy lost in GeV
            new_energy = (energy - energy_lost) #new energy
            if new_energy < 0.0: #check to make sure new energy isnt negative
                muon[0] = 0; muon[2] = 0; continue #set energy and flux to zero manually
            muon[0] = new_energy #new energy
            muon[3] = muon[3] + 10 #add 10 meters
    return final_muons
