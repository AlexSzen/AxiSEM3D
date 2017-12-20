
import numpy as np

def spz2xyz(s, z, phi) :
    r = np.sqrt(np.power(s,2.)+np.power(z,2.))
    if (r<1e-6) :
         theta = 0.
    else :
        theta = np.arccos(z/r)

    x = r * np.sin((theta)) * np.cos(np.deg2rad(phi))
    y = r * np.sin((theta)) * np.sin(np.deg2rad(phi))
    z = r * np.cos((theta))

    return x, y, z


class WavefieldComputer:

    def __init__(self, wavefield, nu, s, z, sem_mesh):
        self.wvf = wavefield
        self.nu = nu
        self.s = s
        self.z = z
        self.sem_mesh = sem_mesh
        self.nu_sum = np.cumsum(np.concatenate(([0],nu)))[:-1]
        self.r = np.sqrt(np.power(s,2.)+np.power(z,2.))
        self.num_steps = len(wavefield)

    def compute_slice(self, r_inner, r_outer, phi):

        r_outer_closest = self.r[np.abs(self.r-r_outer).argmin()]
        r_inner_closest = self.r[np.abs(self.r-r_inner).argmin()]

        tags_slice = np.argwhere((self.r > r_inner_closest) & (self.r < r_outer_closest) )
        n_points_slice = len(tags_slice)

        wvf_slice = np.zeros((self.num_steps,n_points_slice), dtype = np.float32)
        x_slice = np.zeros(n_points_slice, dtype = np.float32)
        y_slice = np.zeros(n_points_slice, dtype = np.float32)
        z_slice = np.zeros(n_points_slice, dtype = np.float32)

        for ip in range(n_points_slice):
            tag = tags_slice[ip]
            elem_tag = np.argwhere(self.sem_mesh==tag)[0] #gives elem, ipol and jpol of this tag. need to solve duplicate points issue
            nuelem = self.nu[elem_tag[0]] #nu of this elem
            ind = self.nu_sum[elem_tag[0]] #index in elem and fouriers dimension of the wavefield
            x_slice[ip], y_slice[ip], z_slice[ip] = spz2xyz(self.s[tag],self.z[tag], phi)

            wvf_slice[:, ip] = self.wvf[:,ind,0,elem_tag[1],elem_tag[2]]
            for ialpha in range(1,nuelem):
                expval = 2 * np.exp(np.deg2rad(phi)*ialpha*1j)
                wvf_slice[:, ip] += np.real(expval * (self.wvf[:,ind+ialpha,0,elem_tag[1],elem_tag[2]] + 1j * self.wvf[:,ind+ialpha,1,elem_tag[1],elem_tag[2]]))

        return x_slice, y_slice, z_slice, wvf_slice

    def compute_shell(self, r_shell, r_tolerance, sample_density, phis, inner_outer):

            tags_shell= np.argwhere(np.abs(self.r-r_shell) < r_tolerance )
            n_arc_points = len(tags_shell)
            if (inner_outer == 'inner'):
                perimeters = ((phis[1]-phis[0])/360.) * 2 * np.pi * self.s[tags_shell] * 1e-6 #in thousand km
            elif (inner_outer == 'outer'):
                perimeters = ((360.-(phis[1]-phis[0]))/360.) * 2 * np.pi * self.s[tags_shell] * 1e-6 #in thousand kilometers
            else :
                raise ValueError("Compute shell : unknown type of shell requested")

            n_shell_points = (perimeters*sample_density).astype(int).squeeze()
            write_offsets = np.cumsum(np.concatenate(([0],n_shell_points)))[:-1]
            tot_shell_points = np.sum(n_shell_points)


            wvf_shell = np.zeros((self.num_steps,tot_shell_points), dtype = np.float32)
            x_shell = np.zeros(tot_shell_points)
            y_shell = np.zeros(tot_shell_points)
            z_shell = np.zeros(tot_shell_points)


            for ip in range(n_arc_points):
                tag = tags_shell[ip]
                elem_tag = np.argwhere(self.sem_mesh==tag)[0]
                nuelem = self.nu[elem_tag[0]]
                ind = self.nu_sum[elem_tag[0]]

                for isamp in range(n_shell_points[ip]):
                    if (inner_outer == 'inner'):
                        phi = phis[0] + isamp * (phis[1]-phis[0]) / n_shell_points[ip]
                    elif (inner_outer == 'outer'):
                        phi = phis[1] + isamp * (360.-(phis[1]-phis[0])) / n_shell_points[ip]

                    x_shell[write_offsets[ip]+isamp], y_shell[write_offsets[ip]+isamp], z_shell[write_offsets[ip]+isamp] = spz2xyz(self.s[tag], self.z[tag], phi)

                    wvf_shell[:,write_offsets[ip]+isamp] = self.wvf[:,ind,0,elem_tag[1],elem_tag[2]]
                    for ialpha in range(1,nuelem):
                        expval = 2 * np.exp(np.deg2rad(phi)*ialpha*1j)
                        wvf_shell[:, write_offsets[ip]+isamp] += np.real(expval * (self.wvf[:,ind+ialpha,0,elem_tag[1],elem_tag[2]] + 1j * self.wvf[:,ind+ialpha,1,elem_tag[1],elem_tag[2]]))

            return x_shell, y_shell, z_shell, wvf_shell
