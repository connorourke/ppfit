import numpy as np
from ppfit.io import read_from_file

def read_optimisation_options( options_file ):
### Read parameters for the optimization routines
    opt_file = read_from_file( options_file, 2, True )
    opt_options = {}
    for i in opt_file:
        if str(i[0]) in ('method_Min','method_BH'):
            opt_options[i[0]] = i[1]
        elif i[0] == 'disp':
            opt_options[i[0]] = bool(i[1])
        elif i[0] in ('maxiter_Min','maxiter_BH','niter_success','niter_BH'):
            opt_options[i[0]] = np.int(i[1])
        else:
            opt_options[i[0]] = np.float(i[1])
    return opt_options

