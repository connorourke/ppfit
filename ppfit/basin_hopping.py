import numpy as np

# These are the bounds for BH
class MyBounds:

    def __init__(self, xmax, xmin ):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        if not tmax and tmin:
            outfile = open('OUTPUT','a')
            outfile.write('Values out of bounds!\n')
            outfile.close()
        return tmax and tmin

# Routine used to write a restart file for accepted parameter sets from BH
# TODO: This should probably be included in the Parameter_Set class, as various output methods that can be called using `callback =` during minimisation
class WriteRestart:

    def __init__( self, nvars, nconst, zerosAndones, nmin, nmax, steps, filename ):
        self.nvars = nvars
        self.nconst = nconst
        self.zerosAndones = zerosAndones
        self.nmax = nmax
        self.nmin = nmin
        self.steps = steps
        self.filename = str( filename )

    def write_bh_restart( self, x, f, accepted ):
        if int( accepted ) == 1: # What's wrong with `if accepted:` ? Is accepted a boolean flag?
            x = np.array(x)
            with open( self.filename, "w") as fl:
                fmt="{0:.7f}"
                fl.write('# The total chi sq is: '+fmt.format(f)+'\n')
                self.write_vars( fl, x )

    def write_local_restart( self, x ):
        x = np.array(x)
        with open( self.filename, "w") as fl:
            self.write_vars( fl, x )

    def write_vars( self, fl, x ):
        fmt="{0:.7f}"
        for lineNumber,var in enumerate(self.nvars):
            output = str( var ) + '\t' + fmt.format( (np.concatenate((self.nconst,x),axis=0) )[lineNumber])+'\t'+ \
                     str( self.zerosAndones[lineNumber] )+'\t'+fmt.format(self.nmin[lineNumber])+'\t'+fmt.format(self.nmax[lineNumber])+'\t' + \
                     fmt.format(self.steps[lineNumber])+'\n'
            fl.write(output)
    

# this routine defines the magnitude of the steps in BH - these steps are defined in PARAMS
class MyTakeStep:

    def __init__(self, stepsizes):
        self.stepsizes = stepsizes

    def __call__(self, x):
        for i in range(len(self.stepsizes)):
            s = self.stepsizes
            x[i] += np.random.uniform(-s[i], s[i])
        return x
