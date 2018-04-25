'''

inputs : unperturbed model parameters, perturbed model params,
kernels obtained in perturbed model

purpose : compute misfit by integrating kernels*perturbation over domain.

'''





import numpy as np
from netCDF4 import Dataset





##### INPUTS #####

INPUT_WVF = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/wavefields/wavefield_db_fwd.nc4'
INPUT_KER = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/kernels/kernels_db.nc4'


##### LOAD FILES #####

dset_fwd = Dataset(INPUT_WVF)
dset_ker = Dataset(INPUT_KER)

dm = dset_fwd.variables['vp'][:] - dset_fwd.variables['vp1D'][:]
m = dset_fwd.variables['vp1D'][:]
kernels = dset_ker.variables["Kernels"][0]
nu = dset_fwd.variables['Nus'][:]
nuKer = dset_ker.variables['Nus'][:]
nr = dset_fwd.variables['Nrs'][:]
nrKer = dset_ker.variables['Nrs'][:]
int_fact = dset_fwd.variables['integral_factor'][:]
nelem = len(nu)


##### FFT, MULTIPLY, INTEGRATE ####

nulineKer = 0
nulineFwd = 0
chi_tot = 0

for ielem in range(nelem):
    
         
    nuelem = min(nu[ielem],nuKer[ielem])
    nrelem = min(nr[ielem],nrKer[ielem])
    chi_elem = 0
        
    ker_c = kernels[nulineKer:nulineKer+nuelem, 4, :, :] + 1j * kernels[nulineKer:nulineKer+nuelem, 5, :, :]

    dm_c = dm[nulineFwd:nulineFwd + nuelem, 0, :, :] + 1j * dm[nulineFwd:nulineFwd+nuelem, 1, :, :]
    int_elem = 2*np.pi*int_fact[ielem,:,:]
    for inr in range(nuelem):
        chi_elem += np.sum( 2*np.real(np.multiply(int_elem,np.multiply(ker_c[inr,:,:], dm_c[inr,:,:])))) #because negative nu are complex conjugate of positive, it simplifies to 2*np.real()

    chi_tot +=(chi_elem)
    nulineKer += nuKer[ielem]
    nulineFwd += nu[ielem]

print(chi_tot/1e15)

