'''
given a kernel with a fourier expansion, looks at
decay of fourier coefficients to produce a new fourier expansion
sufficient to model the kernel.
objective : learn the maximal necessary expansion for kernels.
'''
import numpy as np

class wisdom:

    def __init__(self, kernel, nuKer, epsilon):

        self.kernel = np.squeeze(kernel)
        self.nuKer = nuKer 
        self.nelem = len(nuKer)
        self.nu_sum = np.cumsum(np.concatenate(([0],nuKer)))[:-1]
        self.epsilon = epsilon

    def learn_wisdom(self,comp):
        
        wisdom_nu = []        
        totalMax = np.max(np.abs(self.kernel[:, comp, :, :] + 1j * self.kernel[:, comp+1, :, :]))
        for ielem in range(self.nelem):

            nuElem = self.nuKer[ielem]
            nuOffset = self.nu_sum[ielem]
            wisdom_nu_elem = []
            for ipol in range(5):
                for jpol in range(5):
                    
                    absKerPoint = np.abs(self.kernel[nuOffset:nuOffset+nuElem, comp, ipol, jpol] + 1j * self.kernel[nuOffset:nuOffset+nuElem, comp+1, ipol, jpol])
                    maxKer = np.amax(absKerPoint)
                    indMaxKer = np.argmax(absKerPoint)
                    if maxKer < 1e-7 * totalMax:
                        wisdom_nu_elem.append(0)
                        continue
                    for inu in range(indMaxKer, nuElem):  
                        if absKerPoint[inu] < self.epsilon * maxKer:
                            wisdom_nu_elem.append(inu)
                            break
            
            if wisdom_nu_elem == []:
                max_wisdom_nu_elem = nuElem
            else :
                max_wisdom_nu_elem = np.max(wisdom_nu_elem)

            wisdom_nu.append(max_wisdom_nu_elem)
            
        return wisdom_nu    


