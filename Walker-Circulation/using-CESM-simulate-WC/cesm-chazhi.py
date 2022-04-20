import Ngl
import Nio
import os
import xarray as xr
import numpy as np
import matplotlib as mpl

# python package for plotting maps, 2D plot, etc.
import psyplot.project as psy

# the next line is only necessary when running within a Jupyter notebook
# and allows to inline plots in the Jupyter notebook
%matplotlib inline

# Set figure size for all our plots
mpl.rcParams['figure.figsize'] = [10., 8.]
# get your username from the environment variable USER
username = os.getenv('USER')
# specify the path where your test simulation is stored
path = '/opt/uio/GEO4962/' + username + '/f2000.T31T31.test/atm/hist/'
filename = path + 'f2000.T31T31.test.cam.h0.0009-01.nc'
print(filename)

#  Open the netCDF file containing the input data.
cfile = Nio.open_file(filename,"r")

#  Define the output pressure levels.
pnew = [850.]

#  Extract the desired variables.
hyam = cfile.variables["hyam"][:]
hybm = cfile.variables["hybm"][:]
T    = (cfile.variables["T"][:,:,:,:])
psrf = (cfile.variables["PS"][:,:,:])
P0mb =  0.01*cfile.variables["P0"].get_value()

lats = cfile.variables["lat"][:]
lons = cfile.variables["lon"][:]

#  Do the interpolation.
intyp = 1                              # 1=linear, 2=log, 3=log-log
kxtrp = False                          # True=extrapolate (when the output pressure level is outside of the range of psrf)
  
Tnew = Ngl.vinth2p(T,hyam,hybm,pnew,psrf,intyp,P0mb,1,kxtrp)

Tnew[Tnew==1e30] = np.NaN

T850=xr.Dataset({'T850': (('lats','lons'), Tnew[0,0,:,:])},
                {'lons':  lons, 'lats':  lats})

# plot using psyplot

psy.plot.mapplot(T850, name='T850', title='Temperature (K) at 850 mb')
