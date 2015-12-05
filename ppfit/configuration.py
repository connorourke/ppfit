import sys
import re
import os
from shutil import copyfile
from glob import glob
from ppfit.fitting_parameter import Fitting_Parameter
from ppfit.fitting_data import Forces_Data, Dipoles_Data, Stresses_Data
import numpy as np

def fitting_params_from_fitabinitioin( filename = 'fitabinitio.in' ):
    '''
    Parses 'fitabinitio.in' to obtain the fitting parameters to be adjusted in the fitting procedure.

    Args:
        filename (string) (default 'fitabinitio.in' ): Filename to read fitting parameters from in `fitabinitio` format

    Returns:
        a list of ppfit.Fitting_Parameter objects.
    '''
    with open( filename, 'r') as f:
        data = f.read()
    fitting_params = []
    for line in re.findall( r"MINUIT : fit to ab initio data\n([\s+\w+\n\.-]*)\nPRINTOUT", data )[0].split("\n"):
        if line:
            string, initial_value, max_delta, min_value, max_value = line.split()[1:6]
            fitting_params.append( Fitting_Parameter( string, float( initial_value ), float( max_delta ), float( min_value ), float( max_value ) ) )
    return fitting_params  

class Configuration:

    def __init__( self, runtime_file, restart_file, forces_file, dipoles_file = None, stresses_file = None, nsupercell = 1 ):
        self.runtime = runtime_file
        self.restart = restart_file
        self.training_data = {}
        self.training_data[ 'forces' ] = Forces_Data.load( forces_file )
        if dipoles_file:
            self.training_data[ 'dipoles' ] = Dipoles_Data.load( dipoles_file )
        if stresses_file:
            self.training_data[ 'stresses' ] = Stresses_Data.load( stresses_file )
        self.nsupercell = nsupercell

    @property
    def reference_forces( self ):
        return self.training_data[ 'forces' ].data

    @property
    def reference_dipoles( self ):
        return self.training_data[ 'dipoles' ].data

    @property
    def reference_stresses( self ):
        return self.training_data[ 'stresses' ].data

    def pimaim_run( self ):
        copyfile( self.runtime, 'runtime.inpt' )
        copyfile( self.restart, 'restart.dat' )
        to_delete = glob( '*out*' ) + glob( '*.fort' )
        for f in to_delete:
            os.remove( f )
        os.system( 'pimaim_serial > out.out' )
        # TODO not all of these will be present, depending on the type of calculation
        # TODO can either set this through a PIM / DIPPIM etc. flag, or check whether the files exist
        self.new_forces = np.loadtxt( 'forces.out' )[0::self.nsupercell]
        number_of_ions = self.new_forces.shape[0]
        self.new_dipoles = np.loadtxt( 'dipoles.out' )[0:number_of_ions:self.nsupercell]
        diag_stresses = np.loadtxt( 'xxyyzzstress.out' )[1:4]
        off_diag_stresses = np.loadtxt( 'xyxzyzstress.out' )[1:4]
        self.new_stresses = np.concatenate( ( diag_stresses, off_diag_stresses ) )
        # TODO How should the stress tensors be treated if we have a supercell
       
    #def append_forces( self, dft_force_filename, md_force_filename ):
    #    with open( dft_force_filename, 'a' ) as f:
    #        f.write( ''.join( self.reference_forces ) )
    #    with open( md_force_filename, 'a' ) as f:
    #        f.write( ''.join( self.new_forces ) )

