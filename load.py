import numpy as np
import h5py
import xarray as xr
import devilpy as dp

N_hor = 64

#------
# Read files
#flow = dp.input.read_flowfiles('.')
run = dp.input.read_input('.')

gridfile = h5py.File('grid.h5', 'r')
movfile = h5py.File('movie.h5', 'r')
mov = dp.input.hdf5tods(movfile)
tkefile = h5py.File('tke.h5', 'r')
#meanfile = h5py.File('mean.h5', 'r')
#------

#----
# Set up grids
grids = gridfile['grids']
X_f = np.linspace(0, run.LX, N_hor+1)
Z_f = np.linspace(0, run.LZ, N_hor+1)
Y_f = grids['y'][:]

X, Z = dp.utils.get_grid_centers(X_f, Z_f)
Y = Y_f[:-1]
#----

mov = mov.assign_coords(X=X, Y=Y, Z=Z, label=range(len(mov.label)))

exit()
label = '0001'
w_xz = create_2dda(movfile['w_xy'][label][:], X_c, Y_c, dims='xz')
w_xy = create_2dda(movfile['w_xz'][label][:], X_c, Z_c, dims='xy')

u_xz = create_2dda(movfile['u_xy'][label][:], X_c, Y_c, dims='xz')
u_xy = create_2dda(movfile['u_xz'][label][:], X_c, Z_c, dims='xy')

v_xz = create_2dda(movfile['v_xy'][label][:], X_c, Y_c, dims='xz')
v_xy = create_2dda(movfile['v_xz'][label][:], X_c, Z_c, dims='xy')

T1_xz = create_2dda(movfile['th1_xy'][label][:], X_c, Y_c, dims='xz')
T1_xy = create_2dda(movfile['th1_xz'][label][:], X_c, Z_c, dims='xy')
