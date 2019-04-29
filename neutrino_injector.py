import numpy as np
import time, glob
import fluxes as flux

import argparse
parser = argparse.ArgumentParser(description = "inject neutrinos into our target medium at certain distances")
parser.add_argument('-nu_type', default="solar",type=str, dest = 'nu_type')
parser.add_argument('-num',     default=1e6,    type=int, dest = 'num')
parser.add_argument('-dist',    default=100,    type=int, dest = 'dist')
args = parser.parse_args()

start_time = time.asctime()
print 'Started:', start_time

#array will contain the following entries
#energy, distance, depth, flux, weight
nu_arrray=flux.simulate_neutrinos(args.nu_type, int(args.num), args.dist)

#save the simulation
np.save("data/%s_dist_%s.npy"%(args.nu_type, int(args.dist) ), nu_arrray)
end_time = time.asctime()
print 'Ends:', end_time
