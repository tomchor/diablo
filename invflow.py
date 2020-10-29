import numpy as np
import h5py as h5
import xarray as xr
import devilpy as dp
from glob import glob

outfiles = np.array(sorted(glob("out.*.h5")))[[0,-1]]

run = dp.input.read_grid()

#------
# Read files
allflow = dp.input.read_flowfiles(outfiles=['out.000000.h5', 'out.001000.h5'], rundir=None, variables=['V',
                                                                             'W',
                                                                             'TH1'])
#------

#----
allflow = allflow.assign_coords(X=run.X, Y=run.Y, Z=run.Z)
flow = allflow.isel(itime=0)
