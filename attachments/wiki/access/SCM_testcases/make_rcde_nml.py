# Created by Rob Warren (rob.warren@monash.edu)

# Creates Fortran 90 namelist for running UM Single Column Model in
# radiative-convective equilibrium (RCE) or radiative-convective-dynamical
# equilibrium (RCDE)

# User must specifiy the output file name, a UM vertical levels namelist, the
# run length and model time step, and four parameters relating to the initial
# temperature, specific humidity, and vertical velocity profiles

import f90nml as nml
from f90nml.namelist import Namelist
import numpy as np

# Output namelist file
nmlfile = 'namelist.scm'

# Vertical levels namelist
levsfile = 'vertlevs_L85_50t_35s_85km'

# User-specified run parameters
nday = 100       # run length (days)
dt   = 1200.0    # model time step (s)

# User-specified sounding parameters
t0   = 300.0     # surface temperature (K)
q0   = 18.65e-3  # surface specific humidity (kg/kg)
wmax = 0.05      # maximum vertical velocity (m/s)
hw   = 12500.0   # depth of vertical velocity profile (m)

# Additional sounding parameters
p0   = 101480.0  # surface pressure (Pa)
zt   = 15000.0   # height of tropopause (m)
zq1  = 4000.0    # height scale 1 for specific humidity profile (m)
zq2  = 7500.0    # height scale 2 for specific humidity profile (m)
lr   = 0.0067    # lapse rate (K/m)
qt   = 1.0e-14   # specific humidity above tropopause (kg/kg)

# Other parameters
lon    = 0.0     # longitude (deg)
lat    = 0.0     # latitude (deg)
albedo = 0.07    # surface albedo
co2    = 348     # co2 concentration (ppmv)

# Physical constants
rd = 287.04      # gas constant for dry air (J/kg/K)
g  = 9.79764     # acceleration due to gravity (m2/s2)
cp = 1004.64     # specific heat capacity for dry air (J/kg/K)

########################################################################

# Read in vertical levels namelist
levs = nml.read(levsfile)
ztop = levs['vertlevs']['z_top_of_model']
eta_theta = np.array(levs['vertlevs']['eta_theta'])
eta_rho = np.array(levs['vertlevs']['eta_rho'])

# Convert eta levels to physical heights
z_theta = ztop*eta_theta
z_rho = ztop*eta_rho

# Note the number of levels
nz_theta = len(z_theta)
nz_rho = len(z_rho)
nlev = nz_rho

# Get boolean indices of points in troposphere and stratosphere
itr_theta = (z_theta < zt)
ist_theta = (z_theta >= zt)
itr_rho = (z_rho < zt)
ist_rho = (z_rho >= zt)

# Create specific humidity profile (theta levels)
q = np.zeros(nz_theta)
q[itr_theta] = q0*np.exp(-1*z_theta[itr_theta]/zq1)*np.exp(-1*(z_theta[itr_theta]/zq2)**2)
q[ist_theta] = qt

# Compute virtual temperature at surface and tropopause
tv0 = t0*(1+0.608*q0)
tvt = tv0-lr*zt

# Create virtual temperature profile (theta levels)
tv = np.zeros(nz_theta)
tv[itr_theta] = tv0-lr*z_theta[itr_theta]
tv[ist_theta] = tvt

# Create temperature profile (theta levels)
t = tv/(1+0.608*q)

# Compute pressure at tropopause
pt = p0*(tvt/tv0)**(g/(rd*lr))

# Create pressure profile on theta levels
p_theta = np.zeros(nz_theta)
p_theta[itr_theta] = p0*((tv0-lr*z_theta[itr_theta])/tv0)**(g/(rd*lr))
p_theta[ist_theta] = pt*np.exp(-1*(g*(z_theta[ist_theta]-zt)/(rd*tvt)))

# Compute potential temperature profile
th = t*(100000.0/p_theta)**(rd/cp)

# Create pressure profile (rho levels)
p = np.zeros(nz_rho)
p[itr_rho] = p0*((tv0-lr*z_rho[itr_rho])/tv0)**(g/(rd*lr))
p[ist_rho] = pt*np.exp(-1*(g*(z_rho[ist_rho]-zt)/(rd*tvt)))

# Create vertical velocity profile (theta levels)
w = wmax*np.sin(np.pi*(z_theta)/hw)
w[z_theta == 0] = 0.0
w[z_theta >= hw] = 0.0

# Create ozone profile
ozone = (1e-6)*3.6478*((p/100)**0.83209)*np.exp(-1*(p/100)/11.3515)   # ppmv
ozone *= (48.00/28.96)   # kg/kg

# Convert CO2 concentration to kg/kg
co2 *= (1e-6)*(44.01/28.96)

# Create wind profiles
u = np.zeros(nz_rho)
v = np.zeros(nz_rho)

# Add extra value at top of pressure profile
p = np.append(p, 0.0)

# Set number of forcing times and forcing period
nfor = 2
obs_pd = (nday*86400.)/(nfor-1)

# Create namelist
new = Namelist({})

cats = ['cntlscm',
        'indata',
        'rundata',
        'logic',
        'ingwd',
        'injules',
        'inprof',
        'ingeofor',
        'physwitch',
        'radcloud',
        'inobsfor']

for cat in cats:
#for ii in range(len(cats)):
#    cat = cats[ii]
    new[cat] = Namelist({})

new['cntlscm']['model_levels_nml'] = nz_rho
new['cntlscm']['land_points']      = 0
new['cntlscm']['nfor']             = nfor

new['indata']['lat']  = lon
new['indata']['long'] = lat

new['rundata']['ndayin']   = nday
new['rundata']['nminin']   = 0
new['rundata']['nsecin']   = 0
new['rundata']['timestep'] = dt
new['rundata']['sice_alb'] = albedo
new['rundata']['co2start'] = co2
new['rundata']['co2end']   = co2
new['rundata']['ozone']    = list(ozone)

new['logic']['obs']           = True
new['logic']['obs_surf']      = True
new['logic']['land_sea_mask'] = False
new['logic']['land_ice_mask'] = False
new['logic']['soil_mask']     = False
new['logic']['local_time']    = False
new['logic']['l_qpos_for']    = True

new['inprof']['tstari'] = t0
new['inprof']['kill_interp'] = True   # Model will crash without this!
new['inprof']['z_tom_nml'] = ztop
new['inprof']['p_in']   = list(p)
new['inprof']['theta']  = list(th[1:nlev+1])   # exclude surface
new['inprof']['qi']     = list(q[1:nlev+1])    # exclude surface
new['inprof']['ui']     = list(u)
new['inprof']['vi']     = list(v)
new['inprof']['wi']     = list(w)
new['inprof']['w_advi'] = list(w)

new['inobsfor']['l_vertadv'] = True
new['inobsfor']['obs_pd'] = obs_pd
new['inobsfor']['obs_bot'] = 0.0
new['inobsfor']['obs_top'] = ztop
new['inobsfor']['rlx_t'] = 0 
new['inobsfor']['rlx_q'] = 0
new['inobsfor']['rlx_u'] = 0
new['inobsfor']['rlx_v'] = 0
new['inobsfor']['rlx_w'] = 3
new['inobsfor']['tstar_forcing'] = list(np.zeros(nfor)+t0)
new['inobsfor']['t_inc']  = list(np.zeros(nfor*nlev))
new['inobsfor']['q_star'] = list(np.zeros(nfor*nlev))
new['inobsfor']['u_inc']  = list(np.zeros(nfor*nlev))
new['inobsfor']['v_inc']  = list(np.zeros(nfor*nlev))
new['inobsfor']['w_inc']  = list(np.zeros(nfor*nlev))

# Write namelist to file
new.write(nmlfile,force=True)
