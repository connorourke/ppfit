import numpy as np

class Fitting_Data:

    def __init__( self, filename, scaling = 1.0 ):
        self.filename = filename
        self.scaling = scaling

class Forces_Data( Fitting_Data ):

    def __init__( self, filename, scaling = 1.0 ):
        self.filename = filename
        self.scaling = scaling
        self.data = np.loadtxt( filename )

class Dipoles_Data( Fitting_Data ):
    
    def __init__( self, filename, scaling = 1.0 ):
        self.filename = filename
        self.scaling = scaling
        self.data = np.loadtxt( filename )[:,1:4]

class Stresses_Data( Fitting_Data ):

    def __init__( self, filename, scaling = 1.0 ):
        self.filename = filename
        self.scaling = scaling
        self.data = np.loadtxt( filename ).reshape( [-1,6] )
