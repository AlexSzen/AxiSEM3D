''' Input for wavefield animation '''

#---------- <GENERAL> ----------
INPUT_FILE = '/home/alex/Desktop/VariousAxiSEM3DVersions/kuangdai_axisem/AxiSEM3D/build/output/wavefields/wavefield_db.nc4'
OUTPUT_DIR = '/home/alex/PhD/vtk/'

START_TSTEP = 0 # Read the wavefield from this time step
END_TSTEP = 100  # To this time step
INT_TSTEP = 1   # With this interval
#---------- </GENERAL> ----------

#---------- <PARAMETERS FOR SLICES> ----------
SLICES = 2 # Number of slices to make from the list below.

PHIS_SLICES = [0., 150.] # Azimuth of the slice.

RMIN = [3.5e6, 3.5e6] # Lower radius of the slice.

RMAX = [7.e6, 7.e6] # Upper radius of the slice. If RMAX>R_EARTH --> RMAX=R_EARTH.

#---------- </PARAMETERS FOR SLICES> ----------

#---------- <PARAMETERS FOR SHELLS> ----------

SHELLS = 0 # Number of shells to make from the list below.

INNER_OUTER = ["inner", "outer"] #Inner or outer shell.

PHIS_SHELLS = [[0., 150.], [0., 150.]] #Start and end azimuths of the shell.

R = [3.5e6, 6.371e6] # Mean radius at which shell is computed. GLL points don't have exact radius, so plotted are all point within R_TOLERANCE of R.

R_TOLERANCE = [100000, 10000.] # Spread of radius. 

SAMPLE_DENSITY = [10, 10] # Number of samples every thousand kilometers.
#---------- </PARAMETERS FOR SHELLS> ----------
