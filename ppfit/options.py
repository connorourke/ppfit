import yaml
import re

class Options:
    def __init__( self, options ):
        self.code = options[ 'calculation' ][ 'code' ]
        self.executable = options [ 'calculation' ][ 'exec' ]
        self.code_mpi = options [ 'calculation' ][ 'code_mpi']
        self.mpi_exec = options [ 'calculation' ][ 'mpi_exec']
        self.mpi_np = options[ 'calculation' ][ 'mpi_np' ]
        self.run_configs = options[ 'calculation' ][ 'run_configs' ]
        self.exec_proc = options[ 'calculation' ][ 'exec_proc' ]
        self.mpi_opts = options[ 'calculation' ][ 'mpi_options' ]
        if self.run_configs == 'mpi4py':
           self.cpu_node = options[ 'calculation' ][ 'no_cpu_node']
           self.no_nodes = options[ 'calculation' ][ 'no_nodes' ]
           self.no_workers = options['calculation' ][ 'no_slaves']

def read_options( filename ):
    """Reads in the options file in YAML format

    Args:
        filename (str): filename of the options file

    Returns:
        config  (dict): potential fitting options, as a dictionary

    Notes:
        fixes the YAML resolver for scientific notation: http://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number    
    """
    loader = yaml.SafeLoader
    loader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    options = yaml.load( open( filename ), Loader = loader )
    return options    

