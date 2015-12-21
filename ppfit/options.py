import yaml

def read_options( filename ):
    config = yaml.safe_load( open( filename ) )
    return config    

