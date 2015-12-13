import sys
import numpy as np
from ppfit.io import read_from_file, output
from ppfit.basin_hopping import MyTakeStep, WriteRestart, MyBounds
from scipy.optimize import basinhopping, minimize

class LBFGSB_Minimizer:

    def __init__( self, opts ):
        self.options = { 'ftol': opts.ftol_min,
                         'gtol': opts.gtol_min,
                         'disp': opts.verbose,
                         'maxiter': opts.maxiter_min,
                         'eps': opts.stepsize_min }

    def minimize( self, function, initial_values, bounds ): # can initial values and bounds be passed in as a Fitting_Parameter_Set object?
        print( 'L-BFGS-B minimisation' )
        results_min = minimize( function, initial_values, method = 'L-BFGS-B', bounds = bounds, options = self.options )
        return results_min 

class Nelder_Mead_Minimizer:

    def __init__( self, opts ):
        self.options= { 'ftol': opts.ftol_min,
                        'xtol': opts.xtol_min,
                        'disp': opts.verbose,
                        'maxiter': opts.maxiter_min }

    def minimize( self, function, initial_values ):
        print( 'Nelder-Mead minimisation' )
        results_min = minimize( function, initial_values, method = 'Nelder-Mead', options = self.options )
        return results_min

class CG_Minimizer:

    def __init__( self, opts ):
        self.options = { 'gtol': opts.gtol_min,
                         'disp': opts.verbose,
                         'maxiter': opts.maxiter_min,
                         'eps': opts.stepsize_min }
        self.tol = opts.ftol_min

    def minimize( self, function, initial_values ): # bounds?
        print( 'CG minimisation' )
        results_min = minimize( function, initial_values, method = 'CG', tol = self.tol, options = self.options )
        return results_min
 
def optimise( function, fitting_parameters, opts ):
    method = opts[ 'method' ]
    tot_vars = ( fitting_parameters.to_fit + fitting_parameters.fixed ).strings
    pot_vars = fitting_parameters.to_fit.strings
    const_vars = fitting_parameters.fixed.strings
    const_values = np.asarray( fitting_parameters.fixed.initial_values )
    tot_values_min = ( fitting_parameters.fixed + fitting_parameters.to_fit ).min_bounds
    tot_values_max = ( fitting_parameters.fixed + fitting_parameters.to_fit ).max_bounds
    all_step_sizes = ( fitting_parameters.fixed + fitting_parameters.to_fit ).max_deltas
    step_sizes = np.asarray( fitting_parameters.to_fit.max_deltas )
    to_fit_and_not = ( fitting_parameters.fixed + fitting_parameters.to_fit ).is_fixed
    pot_values_min = fitting_parameters.to_fit.min_bounds
    pot_values_max = fitting_parameters.to_fit.max_bounds
    pot_values = np.asarray( fitting_parameters.to_fit.initial_values )
    minim_bounds = fitting_parameters.to_fit.bounds

    # Choose the calculation order
    if opts[ 'basin_hopping' ][ 'calc_order' ] == 0: 
    # what happens if we are not using basin hopping?
        if opts[ 'method' ] == 'L-BFGS-B':
            minimizer = LBFGSB_Minimizer( opts )
            results_min = minimizer.minimize( function, pot_values, bounds )
        elif opts[ 'method' ] == 'CG':
            minimizer = CG_Minimizer( opts )
            results_min = minimizer.minimize( function, pot_values )
        elif opts[ 'method' ] == 'Nelder-Mead':
            minimizer = Nelder_Mead_Minimizer( opts )
            results_min = minimizer.minimize( function, pot_values )
        else:
            sys.exit( 'minimization method {} not supported'.format( opts[ 'method' ] ) )
        output( results_min.message )
        tot_values = np.concatenate((const_values,results_min.x),axis=0)

        # Write a results file 
        write_Results_min = WriteRestart(tot_vars,const_values,to_fit_and_not,tot_values_min,tot_values_max,all_step_sizes,'RESULTS_min')
        write_Results_min(results_min.x,results_min.fun,accepted=1)

        secondSumOfChi = sumOfChi( potential_file, training_set, chi_squared_scaling, plot = True )
        secondSumOfChi.evaluate(results_min.x)
        mkdir_p('./min-errors-pdfs')
        os.system('mv *.pdf ./min-errors-pdfs')

########################################################################################
# basin hopping part using the minimization parameters as a starting guess             #
# The temperature should be a fraction of the final function value from the minimizer  #
########################################################################################
# Define variables for BH RESTART file 
    write_restart = WriteRestart(tot_vars,const_values,to_fit_and_not,tot_values_min,tot_values_max,all_step_sizes,'RESTART')
#
    if opts[ 'basin_hopping' ]['calc_order' ] == 0:
    # Pass the optimized values to BH
        pot_values = results_min.x
    # Temperature parameter for BH
        temperature = results_min.fun * opts[ 'temperature' ]
    elif opts[ 'basin_hopping' ]['calc_order' ] == 1:
    # Temperature parameter for BH
        temperature = function(pot_values) * opts[ 'basin_hopping' ][ 'temperature' ]
    else:
        exit( 'not recognised as basin hopping calculation order: {}'.format( opts[ 'basin_hopping' ][ 'calc_order' ] ) )
        # Step sizes for BH
    mysteps = MyTakeStep(step_sizes) 
    output( 'The temperature is set to: '+str(temperature)+'\n' )
# Set the options for the minimization algo in BH
    if opts[ 'basin_hopping' ][ 'method' ] in ('L-BFGS-B','CG'):
        options = { 'ftol': opts[ 'basin_hopping' ][ 'tolerance' ][ 'ftol' ],
                    'gtol': opts[ 'basin_hopping' [ 'tolerance' ]][ 'gtol' ],
                    'disp': opts[ 'verbose' ],
                    'maxiter': opts[ 'basin_hopping' ][ 'maxiter' ],
                    'eps': opts[ 'basin_hopping' ][ 'stepsize' ] }
    elif opts[ 'basin_hopping' ][ 'method' ] == 'Nelder-Mead':
        options = { 'ftol': opts[ 'basin_hopping' ][ 'tolerance' ][ 'ftol' ],
                    'xtol': opts[ 'basin_hopping' ][ 'tolerance' ][ 'xtol' ],
                    'disp': opts[ 'verbose' ],
                    'maxiter': opts[ 'basin_hopping' ][ 'maxiter' ] }
    else:
        sys.exit('Minimization method {} not supported'.format( opts[ 'basin_hopping' ][ 'method' ] ) )
# Bounds for BH
    mybounds = MyBounds(pot_values_max,pot_values_min)
    if opts[ 'basin_hopping' ][ 'method' ] in ('L-BFGS-B'):
        # Bounds for minimization method inside BH
        bounds = minim_bounds
        minimizer_kwargs = { 'method': opts[ 'basin_hopping' ][ 'method' ],
                             'bounds': bounds,
                             'options': options}
    else:
        minimizer_kwargs = { 'method': opts[ 'basin_hopping' ][ 'method' ],
                             'options': options}
    results_BH = basinhopping( function, pot_values,
                               T = temperature, 
                               stepsize = opts[ 'basin_hopping'][ 'timestep' ],
                               take_step = mysteps,
                               minimizer_kwargs = minimizer_kwargs,
                               disp = opts[ 'verbose' ],
                               accept_test = mybounds,
                               niter = opts[ 'basin_hopping' ][ 'maxiter' ], # TODO check this
                               callback = write_restart, 
                               niter_success = opts[ 'basin_hopping' ][ 'niter_success' ] )
    output( results_BH.message )
    tot_values = np.concatenate((const_values,results_BH.x),axis=0)
#
## Write a results file for the BH part 
    write_Results_BH = WriteRestart(tot_vars,const_values,to_fit_and_not,tot_values_min,tot_values_max,all_step_sizes,'RESULTS_BH')
    write_Results_BH( results_BH.x, results_BH.fun, accepted = 1 )

    finalSumOfChi = sumOfChi( potential_file, training_set, chi_squared_scaling, plot = True )
    finalSumOfChi.evaluate(results_BH.x)
    mkdir_p('./BH-errors-pdfs')
    os.system('mv *.pdf ./BH-errors-pdfs')
    # close the output file
