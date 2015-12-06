import numpy as np

class Fitting_Data:

    def __init__( self, data ):
        self.data = data

    def __add__( self, other ):
        if kind( self ) is not kind( other ):
            raise( TypeError )
        return Fitting_Data( np.concatenate( self.data, other.data ) )

class Forces_Data( Fitting_Data ):

    @classmethod
    def load( cls, filename ):
        return Forces_Data( data = np.loadtxt( filename ) )
        print( self.id() )
        return self

class Dipoles_Data( Fitting_Data ):
   
    @classmethod 
    def load( cls, filename ):
        return Dipoles_Data( data = np.loadtxt( filename )[:,1:4] )

class Stresses_Data( Fitting_Data ):

    @classmethod
    def load( cls, filename ):
        return Stresses_Data( data = np.loadtxt( filename ).reshape( [-1,6] ) )
