import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

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

class sumOfChi( object ):
  def __init__( self, potential_file, configs ):
    self.potential_file = potential_file
    self.configs = configs

  def __call__( self, test_values ):
    self.potential_file.write_with_parameters( test_values )
   
    for c in self.configs:
        c.pimaim_run()
    
    reference_forces = np.concatenate( [ c.reference_forces for c in self.configs ], axis = 0 )
    new_forces = np.concatenate( [ c.new_forces for c in self.configs ], axis = 0 )

    chiSqF = chi_squared( reference_forces, new_forces, False, 'test' )
    totalChi = chiSqF

    print('totalChi:', totalChi)

    return totalChi
