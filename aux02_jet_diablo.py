import numpy as np
import xarray as xr
from scipy.special import erf

#----
g = 9.81 
α = 2e-4
f0 = 1e-4
π = np.pi
#----

#----
def b_from_W(W, b_inf):
    dbdx = - f0 * W.differentiate('Y')
    b = b_inf + int_indef(dbdx, ivar='X')
    return b
def int_indef(f, ivar):
    Δx = f[ivar].diff(ivar).mean()
    F = f.cumsum(ivar) * Δx
    return F
#----

#-----
X1 = 1600
Y1 = 80
dTdY = 25/1000
#-----

#-----
# From us
X = np.linspace(0, 16e3, 1000)
Y = np.linspace(-500, 500, 1000)

X = xr.DataArray(X, dims='X', coords=dict(X=X))
Y = xr.DataArray(Y, dims='Y', coords=dict(Y=Y))

W0 = -0.85
Y0 = 0
X0 = 8000
T0 = 25/2

f_z = np.exp(-(Y - Y0)**2 / Y1**2)
f_x = erf(-(X-X0)/X1)*np.exp(-(X-X0)**2 / X1**2)

v_cw20 = W0 * f_z * f_x
b_inf = g*(1 + α*dTdY*Y)
b_cw20 = b_from_W(v_cw20, b_inf)

F_x = (1/4)*np.sqrt(π)*X1*(1 - erf(-(X-X0)/X1)**2)
b_anal = b_inf - 2*f0*W0*((Y0-Y)/Y1**2) * f_z * F_x
#-----



#-----
# Absolute vorticity
ζ_cw20 = v_cw20.differentiate('X') + f0
#-----

#-----
N2 = g*α*25/1000
q_cw20 = (ζ_cw20+f0)*N2 - v_cw20.diff('Y')**2*f0
