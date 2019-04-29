import numpy as np
import time, glob
import fluxes as flux

import argparse
parser = argparse.ArgumentParser(description = "inject muons directly into the Earth's surface")
parser.add_argument('-num',     default=1e4,    type=int, dest = 'num')
parser.add_argument('-depth',   default=1e3,    type=int, dest = 'depth')
args = parser.parse_args()

start_time = time.asctime()
print 'Started:', start_time

#energy, flux, weight for muons at surface
muon_arrray_surface=flux.simulate_muons(int(args.num), args.depth)
np.save("data/muons_at_surface.npy",muon_arrray_surface)
muon_arrray_depth=flux.propagate_muons(muon_arrray_surface,args.depth)
np.save("data/muons_at_depth_%s.npy"%args.depth,muon_arrray_surface)
end_time = time.asctime()
print 'Ends:', end_time
