import re

def substitute_parameter( input, to_sub ):
    for k, v in to_sub.items():
        input = re.sub( r"{}".format( k ), str( v ), input )
    return input

class Potential_File:

    def __init__( self, template_filename, potential_parameter_set ):
        with open( template_filename, 'r' ) as file_in:
            self.potential_input = file_in.read()
        self.potential_parameter_set = potential_parameter_set

    # p is a vector of parameter values, that will replace self.potential_parameter_set.to_fit_parameters
    def write_with_parameters( self, p ):
        parameter_values = dict( [ ( k, v ) for k, v in zip( self.potential_parameter_set.to_fit_parameters.strings, p ) ] )
        print( parameter_values )
        potential_to_run = substitute_parameter( self.potential_input, parameter_values )
        with open( 'potential.inpt', 'w' ) as file_out:
           file_out.write( potential_to_run )
