from ppfit.fitting_parameter import Fitting_Parameter

class Fitting_Parameter_Set:

    def __init__( self, fitting_parameters ):
        self.fitting_parameters = fitting_parameters

    @property
    def fixed_parameters( self ):
        return Fitting_Parameter_Set( [ p for p in self.fitting_parameters if p.fixed ] )

    @property
    def to_fit_parameters( self ):
        return Fitting_Parameter_Set( [ p for p in self.fitting_parameters if not p.fixed ] )

    @property
    def strings( self ):
        return [ p.string for p in self.fitting_parameters ]
  
    @property
    def initial_values( self ):
        return [ p.initial_value for p in self.fitting_parameters ]

    @property
    def bounds( self ):
        return [ p.limits for p in self.fitting_parameters ]


    @classmethod
    def from_parameters_file( cls, filename = 'parameters.in' ):
        '''
        Parses 'parameters.in' to obtain the fitting parameters to be adjusted in the fitting procedure.

        Args:
            filename (string) (default 'parameters.in' ): Filename to read fitting parameters from in `parameters` format

        Returns:
            a Fitting_Parameter_Set instance.
        '''
        with open( filename, 'r') as f:
            data = f.readlines()
        fitting_params = []
        for line in data:
            if line[0] != '#':
                string, initial_value, fixed, min_value, max_value, max_delta = line.split()
                fitting_params.append( Fitting_Parameter( string, float( initial_value ), float( max_delta ), float( min_value ), float( max_value ) ) )
        return Fitting_Parameter_Set( fitting_params )
