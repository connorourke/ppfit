import sys
import subprocess
import re
import os
from shutil import copyfile
from glob import glob
from mpi4py import MPI
from ppfit.fitting_parameter import Fitting_Parameter
from ppfit.fitting_data import Forces_Data, Dipoles_Data, Stresses_Data
from ppfit.pimaim_calc import PIMAIM_Run
from socket import gethostname


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

    def __init__( self, options, species, directory, runtime_file, restart_file, forces_file, dipoles_file = None, stresses_file = None, nsupercell = 1 ):
        self.options = options
        self.species = species
        self.parent = os.getcwd()
        self.directory = os.path.join( self.parent, directory )
        self.runtime = runtime_file
        self.restart = restart_file
        self.training_data = {}
        self.training_data[ 'forces' ] = Forces_Data.load( os.path.join( self.directory, forces_file ) )
        if dipoles_file:
            self.training_data[ 'dipoles' ] = Dipoles_Data.load( os.path.join( self.directory, dipoles_file ) )
        if stresses_file:
            self.training_data[ 'stresses' ] = Stresses_Data.load( os.path.join( self.directory, stresses_file ) )
        self.nsupercell = nsupercell

        if self.options.code == 'pimaim':
            self.executable = PIMAIM_Run( self, os.getcwd(), clean = True )
        else:
            sys.exit( '{} not a recognised IP code'.format( code ) )


    @property
    def reference_forces( self ):
        return self.training_data[ 'forces' ].data

    @property
    def reference_dipoles( self ):
        return self.training_data[ 'dipoles' ].data

    @property
    def reference_stresses( self ):
        return self.training_data[ 'stresses' ].data

    @classmethod
    def from_dict( cls,  options, config):
        '''
           Returns configugration object given config and options dicts

        Args:
           config:  dictionary containing the configuration from configs.yml
           options: dictionary containing the options from options.yml
        Returns:
           configuration object
        ''' 
        return cls(      options = options,
                         species = config["species"],
                         directory = config["directory"],
                         runtime_file = config["runtime_file"],
                         restart_file = config["restart_file"],
                         forces_file  = config["forces_file"],
                         dipoles_file = config["dipoles_file"],
                         stresses_file = config["stresses_file"] )

    def run( self, clean = True):
        self.executable.set_up()
        self.executable.run()
        return self.executable.ran_okay

    def collect_data( self, clean = True ):
        self.executable.collect_data()
#        self.executable.tear_down()
