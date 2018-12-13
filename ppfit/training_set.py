import numpy as np
import pathos.multiprocessing as mp
import os
import sys
import time
import subprocess
from socket import gethostname
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()

class Training_Set:

    def __init__( self, configurations, options ):
        self.options = options
        self.configurations = configurations
        self.forces   = np.concatenate( [ c.reference_forces   for c in configurations ] )
        if self.options.dipoles:
            self.dipoles  = np.concatenate( [ c.reference_dipoles  for c in configurations ] )
        self.stresses = np.concatenate( [ c.reference_stresses for c in configurations ] )
        self.cwd = os.getcwd()

        if self.options.run_configs == 'mpi4py':
           if rank < self.options.no_nodes :
                host = open (os.path.join(os.getcwd(),gethostname()),'w')
                host.write(gethostname()+":"+str(self.options.cpu_node))

    def run_config( self, config):
        ran_okay = config.run( clean = True)

        if not ran_okay:
           return( False )
        return( True )

    def run_multi( self ):
        pool_size = int(mp.cpu_count()/self.options.exec_proc)
        pool = mp.Pool(processes=pool_size,maxtasksperchild=1,)
        pool_outputs = pool.map(self.run_config, self.configurations)
        pool.close()
        pool.join()

        for c in self.configurations:
            c.collect_data()
            os.chdir(self.cwd)

        if "False" in pool_outputs:
             return( False )
        else:
             return( True )

    def run_mpi(self):
        MPI.COMM_WORLD.Barrier()

        if rank == 0:
           
           tasks = ([StopIteration] * (self.options.no_workers)) + list(range(0,len(self.configurations)))
           status = MPI.Status()
           while tasks:
             comm.recv(source=MPI.ANY_SOURCE, status=status)
             data=tasks.pop()
             comm.send(obj=data,dest=status.Get_source())
           ran_okay = "True"            
         
 
        if (rank > 0):
           task_data=[]
           cwd=os.getcwd()
           for task in iter(lambda:comm.sendrecv(task_data,dest=0),StopIteration):
                ran_okay = []
                ran_okay.append( self.run_config(self.configurations[task]) )
                os.chdir(cwd)

        MPI.COMM_WORLD.Barrier()
        if "False" in  ran_okay:
           return( False )
        return( True )


    def run( self ):
        if self.options.run_configs == 'serial':
           for c in self.configurations:
              ran_okay = []
              ran_okay.append(self.run_config(c))
              c.collect_data()

        elif self.options.run_configs == 'pool':
           ran_okay = self.run_multi()

           if not ran_okay:
             return( False )
           else:
              return( True  )

        elif self.options.run_configs =='mpi4py':
           ran_okay = self.run_mpi()
           ran_okay = comm.gather(ran_okay, root=0)

           if rank ==  0:
              for c in self.configurations:
                  c.collect_data()
           MPI.COMM_WORLD.Barrier()
           
        if rank == 0:
           if "False" in  ran_okay:
               return( False )
           return( True )
        else:
           return( True )

    @property
    def new_forces( self ):
        return np.concatenate( [ c.new_forces for c in self.configurations ] )

    @property
    def new_dipoles( self ):
        if self.options.dipoles:
            return np.concatenate( [ c.new_dipoles for c in self.configurations ], axis = 0 )

    @property
    def new_stresses( self ):
        return np.vstack( ( c.new_stresses for c in self.configurations ) )


