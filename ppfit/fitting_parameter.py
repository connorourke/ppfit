class Fitting_Parameter:
    
    def __init__( self, string, initial, max_delta, min_value, max_value ):
        self.string = string
        self.initial_value = initial
        self.max_delta = max_delta
        self.limits = ( min_value, max_value )
        self.fixed = self.max_delta == 0.0
