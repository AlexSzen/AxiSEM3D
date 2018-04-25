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
m = dset_fwd.variables['vp'][:]
kernels = dset_ker.variables["Kernels"][0]
nu = dset_fwd.variables['Nus'][:]
nuKer = dset_ker.variables['Nus'][:]
nelem = len(nu)
##### FFT, MULTIPLY, INTEGRATE ####

nulineKer = 0
nulineFwd = 0
chi_tot = 0

for ielem in range(nelem):
    
     
    nuelem = min(nu[ielem],nuKer[ielem])
    test = np.zeros((5,5), dtype = complex)
    chi_elem = 0
    for inu in range(nuelem):
        dm_c = dm[nulineFwd + inu, 0, :, :] + 1j * dm[nulineFwd + inu, 1, :, :]
        ker_c = kernels[nulineKer+inu, 4, :, :] + 1j * kernels[nulineKer + inu, 5, :, :]
        chi_elem +=np.sum( np.multiply(dm_c, ker_c)  )   
        #test += dm[nulineFwd + inu, 0, :, :] + 1j * dm[nulineFwd + inu, 1, : ,:]    
#        print(np.max(dm[nulineFwd + inu, 0, :, :] + 1j * dm[nulineFwd + inu, 1, :, :]) )
        test += kernels[nulineKer + inu, 4, :, :] + 1j * kernels[nulineKer + inu, 5, : ,:]    
    chi_tot +=(chi_elem)
    nulineKer += nuKer[ielem]
    nulineFwd += nu[ielem]

print(chi_tot/1e15)

