import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

fmt="{0:.7f}"

def chi_squared( ai_vals, ff_vals, genplot, filename ):
    ai_s2 = ai_vals.var( axis = 1 )
    genplot = bool(genplot)
    sqDiff = []
    ai_plot = []
    ff_plot = []
    for i in range(np.shape(ai_vals)[-1]):
        numerator = 0.0
        denominator = 0.0
        ai_val = 0.0
        ff_val = 0.0
        for j,k,l in zip(ai_vals,ff_vals,ai_s2):
            if np.shape(ai_vals)[0] < 4:
                numerator += ((k[i]-j[i])**2)
                denominator += l
            else:
                numerator += ((k[i]-j[i])**2)
                denominator += (j[i])**2
            if genplot:
                ai_val += (j[i])**2
                ff_val += (k[i])**2
        sqDiff.append(numerator/denominator)
        if genplot:
            ai_plot.append(ai_val)
            ff_plot.append(ff_val)
    sqDiff = np.array(sqDiff)
    if genplot:
        ai_plot = np.sqrt(ai_plot)
        ff_plot = np.sqrt(ff_plot)
    chiSq = np.sum(sqDiff)/len(sqDiff)
    if genplot:
        with PdfPages(filename+'.pdf') as pdf:
            f, axarr = plt.subplots(2, sharex=True)
            att1 = {'color': 'black', 'markerfacecolor': None, 'markersize': 2.5,
            'markeredgewidth': 0.5, 'alpha': 1.0, 'marker': 'o',
            'markeredgecolor': 'black','linestyle' : ':'}
            att2 = {'color': 'blue', 'markerfacecolor': None, 'markersize': 2.5,
            'markeredgewidth': 0.5, 'alpha': 1.0, 'marker': 'o',
            'markeredgecolor': 'blue','linestyle' : 'None'}
            axarr[0].plot(ai_plot,**att1)
            axarr[0].plot(ff_plot,**att2)
            axarr[0].set_title('Difference = '+str(chiSq))
            axarr[1].plot(sqDiff,'r')
            plt.ylabel('ERRORS')
            plt.xlabel('index')
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
    return chiSq

# This is the objective function evaluated by the minimization algorithms 
class sumOfChi:
  def __init__( self, potential_file, training_set, scaling, plot = False ):
    self.potential_file = potential_file
    self.scaling = scaling
    self.plot = plot
    self.training_set = training_set
    self.ai_forces = training_set.forces.T
    self.ai_dipoles = training_set.dipoles.T
    self.ai_stresses = training_set.stresses.T

  def evaluate( self, test_values ):
    if type( test_values ) is not np.ndarray:
        raise TypeError
    self.potential_file.write_with_parameters( test_values )
    ran_okay = self.training_set.run()
    if not ran_okay:
        totalChi = 1E10
        outfile.write('Error: likely due to unphysical parameter value\n') 
        outfile.close()
        return totalChi
    ff_forces = self.training_set.new_forces.T
    ff_dipoles = self.training_set.new_dipoles.T
    ff_stresses = self.training_set.new_stresses.T

    chiSq = {}
    chiSq[ 'forces' ]   = chi_squared( self.ai_forces, ff_forces, self.plot, 'forces-errors' )
    chiSq[ 'dipoles' ]  = chi_squared( self.ai_dipoles, ff_dipoles, self.plot, 'dipoles-errors' )
    chiSq[ 'stresses' ] = chi_squared( self.ai_stresses, ff_stresses, self.plot, 'stresses-errors' )
    factorTot = sum( self.scaling.values() )
    totalChi = sum( [ self.scaling[ k ] * chiSq[ k ] for k in chiSq.keys() ] ) / factorTot

    outfile = open('OUTPUT','a')
    outfile.write( 'Forces chi sq: ' + fmt.format(chiSq[ 'forces' ])+'\n')
    outfile.write( 'Dipoles chi sq: ' + fmt.format(chiSq[ 'dipoles' ]) +'\n')
    outfile.write( 'Stresses chi sq: ' + fmt.format(chiSq[ 'stresses' ]) + '\n')
    outfile.write( 'Total chi sq (no factors): ' + fmt.format( np.mean( list( chiSq.values() ) ) ) + '\n' )
    outfile.write('Total chi sq: '+fmt.format(totalChi)+'\n')
    outfile.write('\n')
    outfile.close()

    return totalChi

