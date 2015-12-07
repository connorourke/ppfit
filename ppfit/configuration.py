import sys
import subprocess
import re
import os
from shutil import copyfile
from glob import glob
from ppfit.fitting_parameter import Fitting_Parameter
from ppfit.fitting_data import Forces_Data, Dipoles_Data, Stresses_Data
from ppfit.pimaim_calc import PIMAIM_Run

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

    def __init__( self, species, directory, runtime_file, restart_file, forces_file, dipoles_file = None, stresses_file = None, nsupercell = 1 ):
        self.species = species
        self.directory = directory
        self.runtime = runtime_file 
        self.restart = restart_file
        self.training_data = {}
        self.training_data[ 'forces' ] = Forces_Data.load( os.path.join( self.directory, forces_file ) )
        if dipoles_file:
            self.training_data[ 'dipoles' ] = Dipoles_Data.load( os.path.join( self.directory, dipoles_file ) )
        if stresses_file:
            self.training_data[ 'stresses' ] = Stresses_Data.load( os.path.join( self.directory, stresses_file ) )
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

    def run( self, code = 'pimaim', clean = True ):
        if code == 'pimaim':
            executable = PIMAIM_Run( self, clean = clean )
        else:
            sys.exit( '{} not a recognised IP code'.format( code ) )
        executable.set_up()
        executable.run()
        executable.collect_data()
        executable.tear_down()
        return executable.ran_okay
