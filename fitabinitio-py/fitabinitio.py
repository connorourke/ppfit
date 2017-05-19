#! /usr/bin/env python3

import numpy as np
import os
from mpi4py import MPI

from ppfit.chi import chi_squared, sumOfChi
from ppfit.configuration import Configuration
from ppfit.fitting_parameter_set import Fitting_Parameter_Set
from ppfit.potential_file import Potential_File
from ppfit.training_set import Training_Set
from ppfit.inputoutput import mkdir_p
from ppfit.optimisation import optimise
from ppfit.options import read_options, Options

options_read = read_options( 'options.yml' )

comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()


run_options = Options(options_read)




if rank ==0:
 outfile = open('OUTPUT','w')

 config_read = read_options('configs.yml')
 configs=[ Configuration.from_dict( run_options,config) for config in sorted(config_read.values(),key=lambda x:sorted(x.keys()))]
else:
 configs = None


configs = comm.bcast(configs, root=0)


training_set = Training_Set(configs, run_options )

fitting_parameters = Fitting_Parameter_Set.from_parameters_file( 'PARAMS' )
potential_file = Potential_File( 'template_BaTiO3', fitting_parameters )

chi_squared_scaling = { 'forces':   options_read[ 'scaling' ][ 'forces' ],
                        'dipoles':  options_read[ 'scaling' ][ 'dipoles' ],
                        'stresses': options_read[ 'scaling' ][ 'stresses' ] }

 
sum_of_chi = sumOfChi( potential_file, training_set, chi_squared_scaling )

# # evaluate sum_of_chi with our initial potential parameters, and save a .pdf plot
# sum_of_chi.evaluate( fitting_parameters.to_fit.initial_values, plot = True )
# # TODO move plotting capabilities of sum_of_chi to a .plot method, and allow the target directory to be set here
# mkdir_p('./initial-errors-pdfs')
# os.system('mv *.pdf ./initial-errors-pdfs')

optimise( sum_of_chi.evaluate, fitting_parameters, options_read )

MPI.COMM_WORLD.Barrier()
MPI.Finalize()

