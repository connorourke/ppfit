import numpy as np

class Training_Set:

    def __init__( self, configurations ):
        self.configurations = configurations
        self.forces   = np.concatenate( [ c.reference_forces   for c in configurations ] )
        self.dipoles  = np.concatenate( [ c.reference_dipoles  for c in configurations ] )
        self.stresses = np.concatenate( [ c.reference_stresses for c in configurations ] )
