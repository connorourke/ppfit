import numpy as np
import os
import sys
from shutil import copyfile
from glob import glob
import re
from socket import gethostname


def include_in_potential_file( line, species ):
    if not re.search( '#', line ):
        return True
    else:
        species_decoration = re.compile( '#(.*)' ).findall( line )[0].split()
        if set( species_decoration ).issubset( species ):
            return True
    return False

class PIMAIM_Run:

    def __init__( self, configuration, parent, clean = True ):
        self.configuration = configuration
        self.parent = parent
        self.clean = clean
        self.cwd = os.getcwd()
        self.ran_okay = None

    def set_up( self ):
         
        if self.configuration.options.code_mpi == True:
         self.cmd=self.configuration.options.mpi_exec + ' ' + self.configuration.options.mpi_opts + ' ' + self.parent+ '/' + gethostname() + ' ' + str(self.configuration.options.mpi_np) +' '+ str(self.configuration.options.exec_proc) + ' ' +self.configuration.options.executable
        else:
            self.cmd=self.options.executable

        with open( 'potential.inpt' ) as f:
            lines = f.readlines()
        new_potential_file = os.path.join( self.configuration.directory, 'potential.inpt' )
        with open( new_potential_file, 'w' ) as f:
            [ f.write( l )for l in lines if include_in_potential_file( l, self.configuration.species ) ]
        
        self.clean_dir()
        return self

    def run( self ):
        os.chdir(self.configuration.directory)
        os.system( '{} > OUT.OUT'.format( self.cmd ) )
        cg_error = 'cg failed to converge' in open( 'OUT.OUT' ).read()
        if cg_error:
            self.ran_okay = False
        self.ran_okay = True
        #this needs commenting for mpi
        os.chdir(self.parent)
    
    def clean_dir( self ):
        os.chdir(self.configuration.directory)
        to_delete = glob( '*out*' ) + glob( '*.fort' )
        for f in to_delete:
            os.remove( f )

    def collect_data( self ):
        os.chdir(self.configuration.directory)
        self.configuration.new_forces = np.loadtxt( 'forces.out' )[0::self.configuration.nsupercell]
        number_of_ions = self.configuration.new_forces.shape[0]
        self.configuration.new_dipoles = np.loadtxt( 'dipoles.out' )[0:number_of_ions:self.configuration.nsupercell]

        diag_stresses = np.loadtxt( 'xxyyzzstress.out' )[1:4]
        off_diag_stresses = np.loadtxt( 'xyxzyzstress.out' )[1:4]
        self.configuration.new_stresses = np.concatenate( ( diag_stresses, off_diag_stresses ) )
        os.chdir(self.configuration.parent)
         
        # TODO How should the stress tensors be treated if we have a supercell

    def tear_down( self ):
        if self.clean:
            self.clean_dir()
        

