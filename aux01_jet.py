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
x = np.linspace(-6e3, 6e3, 1000)
z = np.linspace(-1e3, 0, 1000)

x = xr.DataArray(x, dims='x', coords=dict(x=x))
z = xr.DataArray(z, dims='z', coords=dict(z=z))
#----

#----
def T_from_v(v, T_west, vert_coor='z', hor_coor='x'):
    dTdx = f0 * v.differentiate(vert_coor) / (g*α)
    Δx = v[hor_coor].diff(hor_coor).mean()
    T = dTdx.cumsum(hor_coor)*Δx + T_west
    return T
def int_indef(f, ivar):
    Δx = v[ivar].diff('x').mean()
    F = f.cumsum(ivar) * Δx
    return F
#----

#-----
X1 = 1600
Z1 = 80
dTdz = 25/1000
#-----

#-----
# From Jiao and Dewar
v0 = 0.35
Z0 = -500
X0 = 4000
T0 = 25
T_west = T0 + dTdz*z
z_arg = -(z - Z0)**2 / Z1**2

v_left = 0*(x+z).where(x<-X0, 0)
v_mid =   (v0 * ((x+X0)/X0)           * np.exp(z_arg)).where((x>=-X0) & (x<0), 0)
v_right = (v0 * np.exp(-x**2 / X1**2) * np.exp(z_arg)).where(x>=0, 0)

v_jd15 = xr.ufuncs.maximum(xr.ufuncs.maximum(v_left, v_mid), v_right)
T_jd15 = T_from_v(v_jd15, T_west)
#-----

#-----
# From us
X = np.linspace(0, 16e3, 1000)
Z = np.linspace(-500, 500, 1000)

X = xr.DataArray(X, dims='X', coords=dict(X=X))
Z = xr.DataArray(Z, dims='Z', coords=dict(Z=Z))

W0 = 0.85
Z0 = 0
X0 = 8000
T0 = 25/2
T_west = T0 + dTdz*Z

f_z = np.exp(-(Z - Z0)**2 / Z1**2)
f_x = erf(-(X-X0)/X1)*np.exp(-(X-X0)**2 / X1**2)

v_cw20 = W0 * f_z * f_x
T_cw20 = T_from_v(v_cw20, T_west, vert_coor='Z', hor_coor='X')
b_cw20 = -g*(1 - α*T_cw20)

fx = f_x
F = (1/4)*np.sqrt(π)*X1*(1 - erf(-(X-X0)/X1)**2)
b_west = -g*(1-α*T_west)
b_inf = g*(1 + α*dTdz*Z)
b_anal = b_inf + 2*f0*v0*((Z0-Z)/Z1**2) * f_z * F
#-----



#-----
# Absolute vorticity
ζ_jd15 = v_jd15.differentiate('x')
ζ_cw20 = v_cw20.differentiate('X')
#-----

#-----
N2 = g*α*25/1000
q_jd15 = (ζ_jd15+f0)*N2 - v_jd15.diff('z')**2*f0
q_cw20 = (ζ_cw20+f0)*N2 - v_cw20.diff('Z')**2*f0
