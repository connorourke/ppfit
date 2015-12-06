import sys
import numpy as np
from ppfit.io import read_from_file
from ppfit.basin_hopping import MyTakeStep, WriteRestart, MyBounds
from scipy.optimize import basinhopping, minimize

def output( msg ):
    outfile = open('OUTPUT','a')
    outfile.write( msg )
    outfile.close()

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
 
class Optimisation_Options:

    def __init__( self, opt_options ):
        self.calc_order =  int(opt_options['calc_order'])
        self.method_min = opt_options['method_min']
        self.ftol_min = opt_options['ftol_min']
        self.gtol_min = opt_options['gtol_min']
        self.xtol_min = opt_options['xtol_min']
        self.maxiter_min = opt_options['maxiter_min']
        self.stepsize_min = opt_options['stepsize_min']
        self.verbose = opt_options['disp']
        self.temperature = opt_options['temperature']
        self.method_BH = opt_options['method_BH']    
        self.ftol_BH = opt_options['ftol_BH']
        self.gtol_BH = opt_options['gtol_BH']
        self.xtol_BH = opt_options['xtol_BH']    
        self.maxiter_BH = opt_options['maxiter_BH']
        self.stepsize_BH = opt_options['stepsize_BH']
        self.timestep = opt_options['timestep']
        self.niter_BH = opt_options['niter_BH']
        self.niter_success = opt_options['niter_success']
        self.scalingF = opt_options['scalingF']
        self.scalingD = opt_options['scalingD']
        self.scalingS = opt_options['scalingS']

def read_optimisation_options( options_file ):
### Read parameters for the optimization routines
    opt_file = read_from_file( options_file, 2, True )
    opt_options = {}
    for i in opt_file:
        if str(i[0]) in ('method_min','method_BH'):
            opt_options[i[0]] = i[1]
        elif i[0] == 'disp':
            opt_options[i[0]] = bool(i[1])
        elif i[0] in ('maxiter_min','maxiter_BH','niter_success','niter_BH'):
            opt_options[i[0]] = np.int(i[1])
        else:
            opt_options[i[0]] = np.float(i[1])
    opts = Optimisation_Options( opt_options )
    return opts

def optimise( function, fitting_parameters, method, opts ):
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
    if opts.calc_order == 0: 
        if method == 'L-BFGS-B':
            minimizer = LBFGSB_Minimizer( opts )
            results_min = minimizer.minimize( function, pot_values, bounds )
        elif method == 'CG':
            minimizer = CG_Minimizer( opts )
            results_min = minimizer.minimize( function, pot_values )
        elif method == 'Nelder-Mead':
            minimizer = Nelder_Mead_Minimizer( opts )
            results_min = minimizer.minimize( function, pot_values )
        else:
            sys.exit( 'minimization method '+method+' not supported' )
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
    if opts.calc_order == 0:
    # Pass the optimized values to BH
        pot_values = results_min.x
    # Temperature parameter for BH
        temperature = results_min.fun * opts.temperature
    elif opts.calc_order == 1:
    # Temperature parameter for BH
        temperature = function(pot_values) * opts.temperature
    # Step sizes for BH
    mysteps = MyTakeStep(step_sizes) 
    output( 'The temperature is set to: '+str(temperature)+'\n' )
# Set the options for the minimization algo in BH
    if opts.method_BH in ('L-BFGS-B','CG'):
        options = { 'ftol': opts.ftol_BH, 
                    'gtol': opts.gtol_BH,
                    'disp': opts.verbose,
                    'maxiter': opts.maxiter_BH,
                    'eps': opts.stepsize_BH }
    elif opts.method_BH == 'Nelder-Mead':
        options = { 'ftol': opts.ftol_BH, 
                    'xtol': opts.xtol_BH,
                    'disp': opts.verbose,
                    'maxiter': opts.maxiter_BH}
    else:
        sys.exit('Minimization method '+method_BH+' not supported')
# Bounds for BH
    mybounds = MyBounds(pot_values_max,pot_values_min)
    if opts.method_BH in ('L-BFGS-B'):
        # Bounds for minimization method inside BH
        bounds = minim_bounds
        minimizer_kwargs = { 'method': opts.method_BH,
                             'bounds': bounds,
                             'options': options}
    else:
        minimizer_kwargs = { 'method': opts.method_BH,
                             'options': options}
    results_BH = basinhopping( function,pot_values,T=temperature,stepsize=opts.timestep,take_step=mysteps,minimizer_kwargs=minimizer_kwargs,disp=opts.verbose,accept_test=mybounds,niter=opts.niter_BH, callback=write_restart, 
                               niter_success = opts.niter_success )
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
