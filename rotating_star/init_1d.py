import StarKillerMicrophysics
import numpy as np

# initialize starkiller
StarKillerMicrophysics.starkiller_initialization_module.starkiller_initialize("probin")
print('initialized')
eos_type_module = StarKillerMicrophysics.Eos_Type_Module()
eos_module = StarKillerMicrophysics.Eos_Module()
network_module = StarKillerMicrophysics.Network()

Gconst = 6.674e-8

class InitModel(object):

    def __init__(self):
        self.nx = 6400 

        self.TOL = 1.e-10
        self.MAX_ITER = 250 
        self.AD_ITER = 2500
        self.MAX_VARNAME_LENGTH = 80

        self.anelastic_cutoff = 9.e4
        self.smallx = 1.e-10

        self.model_file = "15m_500_sec.txt"
        self.mesa = False

        if (self.mesa):
            idxs = ['dens', 'temp', 'pres', 'entr', 'sndspd', 'mass', 'gradr', 'igrav', 'ibrunt', 'iconv_vel', 'ispec']
        else:
            idxs = ['dens', 'temp', 'pres', 'entr', 'enuc', 'spec']

        self.idx = {name: v for (name, v) in zip(idxs, range(1, len(idxs)+1))}

        self.xmin = 0
        self.xmax = 1.75e10
        self.dx = (self.xmax - self.xmin) / self.nx
        self.low_density_cutoff = 1.e-7

        self.temp_fluff_cutoff = 2.e-7 
        self.temp_fluff = 1.e5

        self.varnames_stored = []
        self.nspec = StarKillerMicrophysics.actual_network.nspec
        self.nvar = self.idx['spec'] - 1 + self.nspec

    def read_mesa(self, filename=None):

        if filename is None:
            filename = self.model_file

    def read_file(self, filename=None):

        if filename is None:
            filename = self.model_file

        with open(filename) as f:
            # count number of lines 
            npts_file = sum([1 for line in f])

            # go back to start and read second line in file to get number of variables 
            f.seek(0)
            f.readline()
            l = f.readline() 
            nvars_file = int(l.split(' ')[-1])

            # subtract header rows 
            npts_file -= (nvars_file + 2)

            print(f'{nvars_file} variables found in the initial model file')
            print(f'{npts_file} points found in the initial model file')

            var_idx_map = {}

            # read in the names of the variables 
            for i in range(nvars_file):
                var_name_file = f.readline().strip()

                # create map of file indices to model indices 
                try:
                    var_idx_map[i] = self.idx[var_name_file]
                except KeyError:
                    var_idx_map[i] = network_module.network_species_index(var_name_file.title())
                except:
                    var_idx_map[i] = -1


            base_r = np.zeros(npts_file)
            base_state = np.zeros((npts_file, self.nvar+self.nspec))

            # read in model data 
            for i, line in enumerate(f):
                variables = [float(v) for v in line.split(' ')]

                base_r[i] = variables[0]

                for j in range(nvars_file):
                    if var_idx_map[j] > 0:
                        base_state[i, j] = variables[var_idx_map[j]]


            return npts_file, base_r, base_state

    @staticmethod
    def interpolate(r, npts, model_r, model_var):
        # find location in the coordinate array where we want to interpolate 
        id = npts
        for i in range(npts):
            if model_r[i] >= r:
                id = i 
                break

        if id == 1:

            slope = (model_var[id+1] - model_var[id]) / (model_r[id+1] - model_r[id])
            inter = slope * (r - model_r[id]) + model_var[id]

        else:
            slope = (model_var[id] - model_var[id-1]) / (model_r[id] - model_r[id-1])
            inter = slope * (r - model_r[id]) + model_var[id]

            # safety check to make sure inteprolate lies within the bounding points
            minvar = min(model_var[id], model_var[id-1])
            maxvar = max(model_var[id], model_var[id-1])
            inter = max(inter, minvar)
            inter = min(inter, maxvar)

        return inter

    def init_1d(self):

        eos_state = eos_type_module.eos_t()

        # create coordinates 
        xznl = self.xmin + np.array(range(self.nx)) * self.dx
        xznr = self.xmin + np.array(range(1, self.nx+1)) * self.dx
        xzn_hse = 0.5 * (xznl + xznr)

        # read in the model 
        if self.mesa:
            self.read_mesa()
        else:
            npts_file, base_r, base_state = self.read_file()

        # put the model onto our new uniform grid
        model_hse = np.zeros((self.nx, self.nvar))
        igood = -1 
        for i in range(self.nx):
            for n in range(self.nvar):
                if xzn_hse[i] < base_r[-1]:
                    model_hse[i, n] = self.interpolate(xzn_hse[i], npts_file, base_r, base_state[:,n])

                    igood = i
                else:
                    model_hse[i] = base_state[-1, n]

        # make sure the species fractions sum to 1
        model_hse[:, self.idx['spec']:] /= np.sum(model_hse[:, self.idx['spec']:], axis=1)[:, np.newaxis]

        ###################################################################
        # iterate to find the central density
        ###################################################################

        # because the MESA model likely begins at a larger radius than our first
        # HSE model zone, simple interpolation will not do a good job.  We want to
        # integrate in from the zone that best matches the first MESA model zone,
        # assumming HSE and constant entropy.

        # find the zone in the uniformly gridded model that corresponds to the
        # first zone of the original model
        ibegin = -1 

        for i in range(self.nx):
            if xzn_hse[i] >= base_r[0]:
                ibegin = i 
                break 

        # store the central density. We will iterate until the central density converges 
        central_density = model_hse[0, self.idx['dens']]

        print(f'interpolated central density = {central_density}')

        M_enclosed = np.zeros(self.nx)

        converged_central_density = False

        for iter_dens in range(self.MAX_ITER):

            # compute the enclosed mass 
            M_enclosed[0] = 4/3 * np.pi * self.dx**3 * model_hse[0, self.idx['dens']]

            for i in range(1, ibegin+1):
                M_enclosed[i] = M_enclosed[i-1] + 4/3 * np.pi * (xznr[i] - xznl[i]) * \
                    (xznr[i]**2 + xznl[i]*xznr[i] + xznl[i]**2) * model_hse[i, self.idx['dens']]

            # now start at ibegin and integrate inwards 

            eos_state.t = model_hse[ibegin, self.idx['temp']]
            eos_state.rho = model_hse[ibegin, self.idx['dens']]
            eos_state.xn = model_hse[ibegin, self.idx['spec']:]

            eos_module.eos(eos_type_module.eos_input_rt, eos_state)

            model_hse[ibegin, self.idx['pres']] = eos_state.p

            entropy_want = eos_state.s 

            for i in range(ibegin-1, -1, -1):
                # use previous zone as the initial guess for the temperature and density
                dens_zone = model_hse[i+1, self.idx['dens']]
                temp_zone = model_hse[i+1, self.idx['temp']]
                xn = model_hse[i+1, self.idx['spec']:]

                g_zone = -Gconst * M_enclosed[i] / xznr[i]**2

                # iteration loop 
                ################

                # start off the Newton loop by saying that the zone has not converged 
                converged_HSE = False 

                for iter in range(self.MAX_ITER):

                    p_want = model_hse[i+1, self.idx['pres']] - self.dx * 0.5 * (dens_zone + model_hse[i+1, self.idx['dens']]) * g_zone 

                    # now we have two functions to zero:
                    #   A = p_want - p(rho,T)
                    #   B = entropy_want - s(rho,T)
                    # We use a two dimensional Taylor expansion and find the deltas
                    # for both density and temperature

                    # (t, rho) -> (p, s)

                    eos_state.t = temp_zone
                    eos_state.rho = dens_zone
                    eos_state.xn = xn 

                    eos_module.eos(eos_type_module.eos_input_rt, eos_state)

                    entropy = eos_state.s 
                    pres_zone = eos_state.p

                    dpt = eos_state.dpdt 
                    dpd = eos_state.dpdr 
                    dst = eos_state.dsdt 
                    dsd = eos_state.dsdr 

                    A = p_want - pres_zone
                    B = entropy_want - entropy 

                    dAdT = -dpt 
                    dAdrho = -0.5 * self.dx * g_zone - dpd 
                    dBdT = -dst 
                    dBdrho = -dsd 

                    dtemp = (B - (dBdrho / dAdrho) * A) / ((dBdrho / dAdrho) * dAdT - dBdT)

                    drho = -(A + dAdT * dtemp) / dAdrho 

                    dens_zone = max(0.9 * dens_zone, min(dens_zone + drho, 1.1 * dens_zone))
                    temp_zone = max(0.9 * temp_zone, min(temp_zone + dtemp, 1.1 * temp_zone))

                    if (abs(drho) < self.TOL * dens_zone and abs(dtemp) < self.TOL * temp_zone):
                        converged_HSE = True 
                        break 

                if not converged_HSE:
                    print(f'Error zone {i} did not converge in init_1d')
                    print('integrate up')
                    print(f'dens_zone, temp_zone = {dens_zone}, {temp_zone}')
                    print(f'p_want = {p_want}')
                    print(f'drho = {drho}')

                # call the EOS one more time for this zone and then go on to the next
                # (t, rho) -> (p, s)

                eos_state.t = temp_zone
                eos_state.rho = dens_zone
                eos_state.xn = xn 

                eos_module.eos(eos_type_module.eos_input_rt, eos_state)

                pres_zone = eos_state.p 
                dpd = eos_state.dpdr 

                # update the thermodynamics in this zone 
                model_hse[i, self.idx['dens']] = dens_zone
                model_hse[i, self.idx['temp']] = temp_zone
                model_hse[i, self.idx['pres']] = pres_zone
                model_hse[i, self.idx['entr']] = eos_state.s

            if abs(model_hse[0, self.idx['dens']] - central_density) < self.TOL * central_density:
                converged_central_density = True 

            central_density = model_hse[0, self.idx['dens']]

        if not converged_central_density:
            print(f'ERROR: central density iterations did not converge')

        
        print(f'converged central density = {central_density}')

        ###################################################################
        # compute the full HSE model using our new central density and temperature, 
        # and the temperature structure as dictated by the MESA model 
        ###################################################################

        print('putting the model into HSE on our grid...')

        # compute the enclosed mass
        M_enclosed[0] = 4/3 * np.pi * self.dx**3 * model_hse[0, self.idx['dens']]

        fluff = False 

        for i in range(1, self.nx):

            # use previous zone as initial guess for T, rho 
            dens_zone = model_hse[i-1, self.idx['dens']]
            temp_zone = model_hse[i-1, self.idx['temp']]
            xn = model_hse[i-1, self.idx['spec']:]



if __name__ == "__main__":

    m = InitModel()
    m.init_1d()