''' Input for wavefield animation '''

#---------- <GENERAL> ----------
INPUT_FILE = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/wavefields/wavefield_db_fwd.nc4'
INPUT_FILE2 = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/kernels/kernels_db.nc4'
OUTPUT_DIR = '/home/alex/Desktop/phd/vtk/kernels/'
OUTPUT_DIR_WISDOM = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/input/'

START_TSTEP = 0 # Read the wavefield from this time step
END_TSTEP = 1  # To this time step
INT_TSTEP = 1   # With this interval

PLOT_NU = 0
WISDOM_KERNEL = 0
EPSILON_WISDOM = 1e-9
PLOT_INTFACT = 0
#---------- </GENERAL> ----------

#---------- <PARAMETERS FOR SLICES> ----------
SLICES = 2 # Number of slices to make from the list below.

PHIS_SLICES = [0, 180] # Azimuth of the slice.

RMIN = [0e6, 0e6] # Lower radius of the slice.

RMAX = [7.e6, 7.e6] # Upper radius of the slice. If RMAX>R_EARTH --> RMAX=R_EARTH.

COMP_SLICES = [4, 4] #only even numbers, as need to get comp (real part) and comp+1 (imag)
#---------- </PARAMETERS FOR SLICES> ----------

#---------- <PARAMETERS FOR SHELLS> ----------

SHELLS = 0 # Number of shells to make from the list below.

INNER_OUTER = ["inner", "outer"] #Inner or outer shell.

PHIS_SHELLS = [[0., 1.], [0., 360.]] #Start and end azimuths of the shell.

R = [3.2e6, 6.371e6] # Mean radius at which shell is computed. GLL points don't have exact radius, so plotted are all point within R_TOLERANCE of R.

R_TOLERANCE = [100000, 10000.] # Spread of radius. 

SAMPLE_DENSITY = [10, 10] # Number of samples every thousand kilometers.

COMP_SHELLS = [4,4]
#---------- </PARAMETERS FOR SHELLS> ----------
