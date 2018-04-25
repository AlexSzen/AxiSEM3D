'''
plot netcdf seismograms from AxiSEM3D
'''




from netCDF4 import Dataset
#import animation_utils
import numpy as np
import matplotlib.pyplot as plt
import sys
import os 


INPUT_DIR = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/stations/'
STATION = 'II.AAK.SPZ'

smgr = Dataset(INPUT_DIR + 'axisem3d_synthetics.nc', 'r')


time = smgr.variables['time_points'][:]
station = smgr.variables[STATION][:]


fig = plt.figure()
plt.title(' Displacement : source 90N 0E, rec 20N 0E')
plt.plot(time, station/1e15)
plt.legend(['S','Phi','Z'])
plt.xlabel('Time (s)')
plt.ylabel('Amplitude (m)')
plt.show()
                    

smgr.close()

