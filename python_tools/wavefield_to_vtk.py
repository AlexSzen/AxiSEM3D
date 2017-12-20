'''
Load wavefield netcdf file and converts to vtk for paraview visualization.
Wavefield has dimensions [time][fouriers and elems][6 : real and imag parts of s p z][ipol][jpol]
Author : Alexandre Szenicer
'''

from netCDF4 import Dataset
from scipy.spatial import KDTree
import pyvtk
import numpy as np
import time
import sys
import os
from wavefield_utils import WavefieldComputer
from params_animation import *
#---------- <BEGIN TIMER> ----------
t0 = time.time()
#---------- </BEGIN TIMER> ----------

#---------- <INPUT AND OUTPUT> ----------

if (SLICES):
    if not os.path.exists(OUTPUT_DIR+"slices"):
        os.makedirs(OUTPUT_DIR+"slices")

if (SHELLS):
    if not os.path.exists(OUTPUT_DIR+"shells"):
        os.makedirs(OUTPUT_DIR+"shells")




#---------- </INPUT AND OUTPUT> ----------

#---------- <LOAD VARIABLES> ----------
print("Start")
dset_fwd = Dataset(INPUT_FILE)

wavefield = dset_fwd.variables["displacement_wavefield"][START_TSTEP:END_TSTEP:INT_TSTEP] #domain decomposed wavefield.
sem_mesh = dset_fwd.variables["sem_mesh"][:] #for each point contains global point tag
s_temp = dset_fwd.variables["mesh_S"][:]
z_temp = dset_fwd.variables["mesh_Z"][:]
nu = dset_fwd.variables['Nus'][:]

npoints = np.max(sem_mesh) + 1
nelem = len(nu)
s = np.zeros(npoints, dtype = np.float32)
z = np.zeros(npoints, dtype = np.float32)
for ielem in range(nelem):
    s[sem_mesh[ielem]] = s_temp[ielem]
    z[sem_mesh[ielem]] = z_temp[ielem]


num_steps = int((END_TSTEP - START_TSTEP)/INT_TSTEP)

#---------- </LOAD VARIABLES> ----------

#---------- <CREATE VTK FILE> ----------

wc = WavefieldComputer(wavefield, nu, s, z, sem_mesh)

for i_slice in range(SLICES):
    x_slice, y_slice, z_slice, wvf_slice = wc.compute_slice(RMIN[i_slice], RMAX[i_slice], PHIS_SLICES[i_slice])
    points_slice = list(zip(x_slice,y_slice,z_slice))
    range_slice = range(len(x_slice))
    for it in range(num_steps):
        vtk = pyvtk.VtkData(
                pyvtk.UnstructuredGrid(points_slice,range_slice),
                pyvtk.PointData(pyvtk.Scalars(wvf_slice[it],name='wavefield_slice')),
                'animation'
                )
        vtk.tofile(OUTPUT_DIR + 'slices/' + 'slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
        + str(PHIS_SLICES[i_slice]) + '_' + str(it) + '.vtk')

    ### write info file for each slice
    f = open(OUTPUT_DIR + 'slices/' 'slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
    + str(PHIS_SLICES[i_slice]) + '_INFO.txt','w')
    f.write('########## SLICE INFO ##########\n' )
    f.write('PHI (degrees) ' + str(PHIS_SLICES[i_slice]) + '\n')
    f.write('RMIN (m) ' + str(RMIN[i_slice]) + '\n')
    f.write('RMAX (m)' + str(RMAX[i_slice]) + '\n')
    f.write('NUMBER_POINTS ' + str(len(x_slice)) + '\n')
    f.close()


for i_shell in range(SHELLS):
    x_shell, y_shell, z_shell, wvf_shell = wc.compute_shell(R[i_shell], R_TOLERANCE[i_shell], SAMPLE_DENSITY[i_shell], np.sort(PHIS_SHELLS[i_shell]), INNER_OUTER[i_shell])
    points_shell = list(zip(x_shell,y_shell,z_shell))
    range_shell = range(len(x_shell))
    for it in range(num_steps):
        vtk = pyvtk.VtkData(
                pyvtk.UnstructuredGrid(points_shell,range_shell),
                pyvtk.PointData(pyvtk.Scalars(wvf_shell[it],name='wavefield_'+INNER_OUTER[i_shell])),
                'animation'
                )
        vtk.tofile(OUTPUT_DIR + 'shells/' + INNER_OUTER[i_shell] + '_shell_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '.vtk')

    ### write info file for each shell
    f = open(OUTPUT_DIR + 'shells/' + INNER_OUTER[i_shell] + '_shell_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
    + '_' + str(PHIS_SHELLS[i_shell][1]) + '_INFO.txt', 'w')
    f.write('########## SHELL INFO ##########\n' )
    f.write('INNER/OUTER ' + INNER_OUTER[i_shell] + '\n')
    f.write('PHIS (degrees) ' + str(PHIS_SHELLS[i_shell]) + '\n' )
    f.write('MEAN_RADIUS (m)' + str(R[i_shell]) + '\n')
    f.write('TOLERANCE_RADIUS (m) ' + str(R_TOLERANCE[i_shell]) + '\n')
    f.write('SAMPLE_DENSITY (m^(-6))' + str(SAMPLE_DENSITY[i_shell]) + '\n')
    f.write('NUMBER_POINTS ' + str(len(x_shell)) + '\n')
    f.close()

#---------- </CREATE VTK FILE> ----------

#---------- <END TIMER> ----------
t2 = time.time()
print("Total runtime is " + str(t2-t0) + "s")
#---------- </END TIMER> ----------
