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
