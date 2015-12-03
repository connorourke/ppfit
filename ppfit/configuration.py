import sys
import re
import os
from shutil import copyfile
from glob import glob
from ppfit.fitting_parameter import Fitting_Parameter

def substitute_parameter( input, to_sub ):
    for k, v in to_sub.items():
        input = re.sub( r"{}".format( k ), str( v ), input )
    return input

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
            fitting_params.append( Fitting_Parameter( string, initial_value, max_delta, min_value, max_value ) )
    return fitting_params  

def initialise_potential_file( potential_file ):
    with open( potential_file, 'r' ) as file_in:
        potential_input = file_in.read()
    to_sub = {}
    potential_fitting_parameters = fitting_params_from_fitabinitioin()
    for i, fit_param in enumerate( potential_fitting_parameters ):
        to_sub[ fit_param.string ] = sys.argv[ i+1 ]
    potential_to_run = substitute_parameter( potential_input, to_sub )
    with open( 'potential.inpt', 'w' ) as file_out:
        file_out.write( potential_to_run )

class Configuration:

    def __init__( self, runtime_file, restart_file, forces_file, nsupercell ):
        self.runtime = runtime_file
        self.restart = restart_file
        with open( forces_file, 'r' ) as f:
            self.reference_forces = f.read()
        self.nsupercell = nsupercell

    def pimaim_run( self ):
        copyfile( self.runtime, 'runtime.inpt' )
        copyfile( self.restart, 'restart.dat' )
        to_delete = glob( '*out*' ) + glob( '*.fort' )
        for f in to_delete:
            os.remove( f )
        os.system( 'pimaim_serial > out.out' )
        with open( 'forces.out', 'r' ) as f:
            self.new_forces = f.readlines()[0::self.nsupercell]

    def append_forces( self, dft_force_filename, md_force_filename ):
        with open( dft_force_filename, 'a' ) as f:
            f.write( ''.join( self.reference_forces ) )
        with open( md_force_filename, 'a' ) as f:
            f.write( ''.join( self.new_forces ) )
