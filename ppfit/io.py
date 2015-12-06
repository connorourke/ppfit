import os, errno

def read_from_file(name,cols,split):
    '''
    This routine reads the input files, e.g. TEMPLATES, OPT, AI_*,etc
    '''
    configs = []
    cols = int(cols)
    split = bool(split)
    with open(name, "r") as f:
        for line in f.readlines():
            li = line.strip()
            if not li.startswith("#"):
                if cols == 1:
                    if split == True:
                        configs.append(line.split())
                    else:
                        configs.append([line])
                if cols == 2:
                    configs.append(line.split())
                if cols == 3:
                    configs.append(line.split())
                if cols == 6:
                    configs.append(line.split())
    return configs

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: 
            raise

def output( msg ):
    outfile = open('OUTPUT','a')
    outfile.write( msg )
    outfile.close()
