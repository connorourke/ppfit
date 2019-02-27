# A Guide to Potential Fitting With PPFIT

A guide to fitting [interaction potentials](http://pubs.rsc.org/en/Content/ArticleLanding/2003/FD/b300319c#!divAbstract) to DFT calculations using [VASP](https://www.vasp.at/) and [ppfit](https://github.com/bjmorgan/ppfit) and  the PIMAIM molecular dynamics code.

There are two main steps in generating an interaction potential with this approach:

- **Generating a training set**: Using VASP to perform calculations on a set of relevant structures with which to train the potential.
- **Potential Fitting**: Using this training set, minimise the errors between the VASP training set data and PIMAIM data produced from the potential by varying the potential parameters.

## Generating Training Sets

The three quantities we will fit our potential to are:

<ul style="list-style-type:disc">
  <li>Atomic forces</li>
  <li>Atomic dipoles</li>
  <li>Stresses on the crystal</li>
</ul>

So we wish to perform VASP calculations on a relevant and broad set of structures to obtain these properties. 

The bigger and more varied the set of structures in the training set, the better and more transferrable the resulting potential will be - however fitting a potential for a large and varied training set will also prove more difficult.

### A few tools for setting up the training set structures

As an example lets consider the case where the end goal is to obatain a potential for modelling how a dopant diffuses through a crystal. We will probably want to consider how this diffusion varies with dopant concentration, so our training set should contain a range of structures at different concentrations.

In order to generate this training set [pymatgen](http://pymatgen.org/)  can be very useful:

```python
from pymatgen import Structure
	
#set up structure from POSCAR file
struct = Structure.from_file("POSCAR")
	
#now we have a structure we can change species on a site:
	
struct[0]= "Li"
```	
	
There are endless possibilities on manipulating structures with pymatgen, so check out the documentation and examples. Also in [crystal_torture](https://github.com/connorourke/crystal_torture) there are few simple doping routines which may be useful. For example taking a clean, relaxed crystal structure and generating a range of doped structures at varying concentrations:

```python	
import crystal_torture.pymatgen_doping as pd
from pymatgen import Structure
import copy
	
#set up a structure from unit cell
struct = Structure.from_file("POSCAR_unit.vasp")
#you can label sites, eg:
struct.add_site_property("label",["A"]*8+["B"]*16+["O"]*32)
	
	
#create a 2x2x2 supercell
struct.make_supercell([2,2,2])
	
#dope at various concentrations and output new structure
# Eg Remove 2 Mg on A-sites and insert 1Al + 1 Li in their place
	
for Li_conc in [0.1,0.2,0.3,0.4]:
    doped_structure = copy.deepcopy(supercell)
    doped_structure = pd.dope_structure(doped_structure,conc=Li_conc,species_to_rem="Mg",species_to_insert=["Li","Al"],label_to_remove="A")
	    dope_structure.to(filename="POSCAR_"+ str(Li_conc) + ".vasp"
```	    
	    
Another useful tool is [bsym](https://github.com/bjmorgan/bsym) , which allows you to generate a set of symmetry-inequivalent crystal structures - and thereby making sure that the same structure doesn't occur twice in the training set.

### Generating Thermally Distorted Structures

First, for each structure in your training set, do a standard calculation to relax the crystal lattice and atoms, ensuring all the usual convergence checks etc have been made (If you haven't used VASP before there are plenty of tutorials online. A starting point can be found on the VASP website [here](http://www.vasp.at/vasp-workshop/tutorials/) ).

Next generate a set of thermally distorted structures from our training set of relaxed crystal structures using MD (again plenty online on doing MD with VASP). Example INCAR:

	SYSTEM = Spinel 

	!Start Parameters:
	 NWRITE = 0        (Low-level output information for MD)   
	 ISTART = 1        (Read existing wavefunction) 
	 INIWAV = 1        (Random initial wavefunction)
	 ICORELEVEL = 1    (Print corelevels in OUTCAR))
	
	!Electronic Relaxation:
	 PREC  =  Low (Precision level)
	 LREAL = Auto      (Projection operators: automatic)
	 ROPT = 1E-03 1E-03 1E-03 1E-03 1E-03
	 ALGO  = Very Fast   (Elect. algorithm for MD)

	 ENMAX = 400.00 eV (Plane-wave cutoff) 
	 NELM  = 20        (Max number of SCF steps)   
	 !NELMIN = 4        (Min number of SCF steps) 
	 EDIFF = 5E-05     (SCF convergence) 
	 ISPIN =  1        (Closed shell)
	 GGA   =  PS       (PBE exchange-correlation)
	 ADDGRID = .TRUE.  (Increase grid: helps GGA convergence) 
	 LASPH   = .TRUE.  (Non-spherical elements: PAW d/f convergence)
	
	!Ionic Relaxation:
	 EDIFFG =  -0.010   (Ionic convergence eV/A)
	 NSW    =  1000    (Max ionic steps) 
	 NBLOCK =     20     (Update XDATCAR/DOSCAR every X steps) 
	 IBRION =      0     (Ions: 0-MD, 1-Quasi-New, 2-CG)
	 ISIF   =      2     (Stress/Relaxation: 2-Ions, 3-Shape/Ions/V, 7-Vol)
	 ISYM   =      0     (Symmetry: Use all, 0: none) 
	 LCORR  =      F     (Add non-SCF force correction)
	 ISMEAR =      0     (Gaussian smearing, Metals:1, MP)
	
	!Miscellaneous:
	 LORBIT    =   11     (PAW radii for projected DOS)
	
	!Molecular Dynamics:
	 POTIM  =  2.0  (Timestep fs)
	 MDALGO = 2 (NH thermostat, define SMASS)
	 TEBEG  =     1500  (Start temp K)
	 TEEND  =     1500  (End temp K)
	 SMASS  =     1  (T scaling every NBLOCK stps in NH NVT)

	 MAXMIX = 40

For a good potential we want suitably large forces, dipoles etc so pick a high temperature (without actually melting the crystal....).

Once each structure is equilibrated we can take a snapshot during MD as our thermally distorted training set structure.

### Getting Forces, Dipoles and Stresses

Once you have a set of structures you wish to generate a training set for, it is on to actually getting the forces, dipoles and stresses out for each structure. 

The dipoles are calculated from maximally-located Wannier functions using the [wannier90](http://wannier.org/) code. There are two ways of using the code - either standalone, or as a library which plugs directly into VASP. The library which plugs into VASP is useful, but no longer maintained and as far as I have found **will only work with the intel compilers**.

There is a version of VASP on Balena which has been linked to wannier90. Example run_script for wannier90 linked VASP on Balena:

	#!/bin/bash

	# set the account to be used for the job
	#SBATCH --account=free

	# set name of job
	#SBATCH --job-name=SpSi19156
	#SBATCH --output=SpSi19156.%j.o
	#SBATCH --error=SpSi19156.%j.e

	# set the number of nodes and partition
	#SBATCH --nodes=2
	#SBATCH --ntasks-per-node=16
	#SBATCH --partition=batch-128gb

	# set max wallclock time
	#SBATCH --time=24:00:00

	# Load dependant modules

	module load untested
	module load vasp/intel/5.3.5-wannier

	# run the calc
	mpirun -np $SLURM_NTASKS vasp > output

Or you can compile wannier90 ([download wannier90](http://wannier.org/code/wannier90-2.1.0.tar.gz), extract and follow the instructions in README.install). An example make.sys for Balena is:


	#=====================================================
	# For Linux with intel version 11/12 on 64bit machines
	#=====================================================
	F90 = ifort
	COMMS=mpi
	MPIF90=mpifort
	FCOPTS=-O2 -heap-arrays
	LDOPTS=-O2
	
	#========================================================
	# Intel mkl libraries. Set LIBPATH if not in default path
	#========================================================
	
	LIBDIR = /opt/intel/mkl/lib/intel64
	LIBS   =  -L$(MKLROOT)/lib/intel64/ -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core -liomp5 -lpthread -lmkl_intel_lp64
	
And to link with VASP include the precompiler flag when comilinng VASP:

	-DVASP2WANNIER90

and change the `LIB` variable in the makefile to:

	LIB     = -L../vasp.5.lib -ldmy  \
	     ../vasp.5.lib/linpack_double.o $SRC_DIR//wannier90-1.2/libwannier.a $(SCA) $(LAPACK) $(BLAS)

where `$SRC_DIR` is the path to src folder conatinng wannier90 folder.

### Running wannier90

`wannier90.x`, the stand alone executable of wannier90, compiled with the flags above, will do multithreading, (there is a full mpi parallel version, postw90.x, but it cannot be used to get the wannier functions). So running wannier90 directly with VASP is not a good idea if you are running on more than 1 node. If you do, once the VASP run finishes, you will be left with processes sitting idle while the wannierisation is done on a single node. Instead use the wannier90 linked version of VASP to setup the input files for `wannier90.x`.

First run VASP to get the `wannier90.x` input file by adding:

	! Wannier90 interface
	 LWANNIER90 = .TRUE.
	 LWANNIER90_RUN = .FALSE.
	 LWRITE_MMN_AMN = .FALSE.

to your INCAR file. This will produce a `wannier90.win` input file for `wannier90.x` that will look something like:

	num_wann =  1152  ! set to NBANDS by VASP
	
	begin unit_cell_cart
	    16.0047880    -0.0125190     0.0062410
	    -0.0125330    15.9861200    -0.0601310
	     0.0062910    -0.0602240    15.9443510
	end unit_cell_cart
	
	begin atoms_cart
	Li       2.7352152     6.6092436     7.3561511
	Li       3.1519070     7.2170067    15.1044698
	Li       1.1158266     1.1396616     4.9965117
	Li       2.7689680    14.9356282    15.2954379
	Li       9.2449171     8.8026930     4.7331089
	Li      11.1098388     7.1236397    14.9664858
	Li       8.9249957     0.8921718     4.7658370
	Li      11.1163191    14.9607525    15.0853954
	Li       5.0173436     1.1283086     0.5202051
	Li       4.7543198     1.0457922     9.2539624
	Li       4.9477804     9.0352103     0.7145018
	Li       5.0039462     8.7704106     8.7980751
	Li      13.2497975     0.8232034     1.1766713
	Li      13.1232487     0.5689080     9.1013687
	Li      13.2758036     9.1587768     1.0812255
	Li      13.0412798     9.0387761     8.7883711
	Mg       3.3157509     2.8938116     2.8106749
	Mg       3.0311046     2.6632085    10.9474274
	Mg       3.0138803    10.8505401     3.2145992
	Mg       2.6096248    10.6228537    11.1981178
	Mg      10.8603323     2.9979929     3.1528677
	Mg      11.0753543     2.9764755    11.1282901
	Mg      11.3314649    10.9836656     3.1106120
	Mg      10.9514092    10.7961792    11.3086173
	Mg       6.9774056     6.9419597     2.7976106
	Mg       7.0210682     6.7782793    10.5817932
	Mg       7.1786750    14.8867728     2.7254015
	Mg       6.8630333    14.9263805    10.6265118
	Mg      14.8408487     7.0135163     2.9094471
	Mg      14.9212948     6.7939681    10.8912722
	Mg      14.9877650    14.8361424     2.9717474
	O        1.9687028     6.1144890     2.0621202
	O        2.2512478     5.9128638    10.0028526
	O        1.9952795    13.9613237     2.1996783
	O        2.1095583    13.9797569    10.1533038
	O       10.2285076     5.8772415     2.2146506
	O       10.3947233     6.0620833     9.9992191
	O       10.2706964    13.9680453     1.9651068
	O       10.1457673    13.9659481    10.0098329
	O        1.9455014     1.8467679     1.8936887
	O        1.8415577     1.8939377     9.6931643
	O        1.9451027     9.9711763     1.9215546
	O        1.9148067     9.7192281     9.6484613
	O        9.8648484     2.0661101     1.8773128
	O        9.9377715     1.9561379     9.8286419
	O       10.0836110    10.0784196     1.9124603
	O        9.9246869     9.9245710     9.9080311
	end atoms_cart
	
	mp_grid =     3     3     3
	
	begin kpoints
	      0.000000000000      0.000000000000      0.000000000000
	      0.333333333333     -0.000000000000      0.000000000000
	     -0.000000000000      0.333333333333     -0.000000000000
	      0.333333333333      0.333333333333      0.000000000000
	     -0.333333333333      0.333333333333     -0.000000000000
	      0.000000000000      0.000000000000      0.333333333333
	      0.333333333333      0.000000000000      0.333333333333
	     -0.333333333333      0.000000000000      0.333333333333
	     -0.000000000000      0.333333333333      0.333333333333
	      0.333333333333      0.333333333333      0.333333333333
	     -0.333333333333      0.333333333333      0.333333333333
	      0.000000000000     -0.333333333333      0.333333333333
	      0.333333333333     -0.333333333333      0.333333333333
	     -0.333333333333     -0.333333333333      0.333333333333
	     -0.333333333333      0.000000000000      0.000000000000
	      0.000000000000     -0.333333333333      0.000000000000
	     -0.333333333333     -0.333333333333      0.000000000000
	      0.333333333333     -0.333333333333      0.000000000000
	      0.000000000000     -0.000000000000     -0.333333333333
	     -0.333333333333     -0.000000000000     -0.333333333333
	      0.333333333333      0.000000000000     -0.333333333333
	      0.000000000000     -0.333333333333     -0.333333333333
	     -0.333333333333     -0.333333333333     -0.333333333333
	      0.333333333333     -0.333333333333     -0.333333333333
	     -0.000000000000      0.333333333333     -0.333333333333
	     -0.333333333333      0.333333333333     -0.333333333333
	      0.333333333333      0.333333333333     -0.333333333333
	end kpoints
	

Now in order to generate the overlap matrices wannier needs for wannierisation we need to run again, but to only include occupied bands. To do this grep for  the number of electrons in the OUTCAR file (and divide by two if spin unpolarised):

	grep NELECT OUTCAR | awk '{print $3/2}'

and edit the wannier90.win file to exclude the unoccupied bands:

	num_wann =  1024  ! set to NBANDS by VASP
	
	exclude_bands = 1025-1152
	
	use_bloch_phases = T
	
	write_xyz = T Write *.xyz file
	num_dump_cycles = 200
	num_iter = 100000
	conv_tol = 0.001
	conv_window = 2

The last lines here set the variables for wannierisation needed in the next step (for more info read the [wannier90 documentation](http://www.wannier.org/doc/user_guide.pdf)).

Running VASP again, but setting `LWRITE_MMN_AMN = .TRUE.` will now produce the overlap file wannier90.mmn. 

Now we can run `wannier90.x` directly on a single node. Note the default filename seed used by `wannier90.x` is `wannier` not `wannier90` as output by vasp, so `wannier90.win` will have to be renamed to `wannier.win` or wannier90.x will need to be run as:

	wannier.x wannier90 > wannier.output

where wannier90 is the seedname.

Throughout the `wannier90` run you can check the convergence by:

```
grep CONV output.wannier
```
which will give output similar to:

``` 
+--------------------------------------------------------------------+<-- CONV
 | Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
 +--------------------------------------------------------------------+<-- CONV
      0     0.681E+03     0.0000000000      680.5027841918       5.51  <-- CONV
      1    -0.103E-01     0.7375245734      680.4925116260     780.59  <-- CONV
      2    -0.132E-01     0.5314444713      680.4792848117    1560.82  <-- CONV
```
The `Delta Spread` is the spread of the wannier functions at that stage in the wannierisation, and it is this property you wish to converge.  It is worth keeping an eye on  `RMS Gradient` too - this will give an indicaton of how close you are to getting a converged set of wannier functions out. Occasionally the wannierisation will stop having hit the convergence criteria while there is still a large RMS gradient, this can indicates the wannierisation hasn't properly converged, so tighten the convergence criteria and run again (setting `conv_window` > 2 should also avoid this problem).

When wannier90.x has finished you should have a `wannier90_centres.xyz` file containing the wannier centres.

To generate the dipoles and extract the forces and stresses use this [wannier2dipoles.py](https://github.com/connorourke/ppfit/blob/develop/ppfit/scripts/wannier2dipoles.py) script:

	wannier2dipoles.py -p POSCAR -wf wannier90_centres.xyz -u Bohr 

(You can also reorder the atoms by type. For info on arguments use `wannier2dipoles.py -h`)

After running wannier2dipoles.py you will be left with four files: `forces.out`, `dipoles.out`, `stresses.out` and `restart.dat`. These files contain the target forces, dipoles and stresses as well as the structure file for your member of the training set.

## Potential Fitting

Now you have a suitable training set, it is on to fitting a potential using [ppfit](https://github.com/bjmorgan/ppfit). This guide will walk through fitting a potential using `ppfit` on Balena.

### Setting up ppfit

`ppfit` takes each configuration in your training set, performs a single MD step with the `pimaim_mpi`  MD code and optimises the potential parameters with respect to the errors between the output values and that of the training set  (though you could use another MD code, and indeed generate your DFT training set with a code other than vasp).

`ppfit` can run in three ways:
- **Serially**: performing the MD on each configuration in turn.
- **Using multiprocessing**: performing the MD steps on each configuration concurrently on a single node.
- **Using mpi4py**: performing the MD steps on each configuration concurrently on multiple nodes.

For the purpose of this guide we will walk through perfroming a fit the quickest way, using `mpi4py` on the University of Bath supercomputer, Balena (the other two approaches are simpler and should follow easily).

First we need a working version of `pimaim_mpi` on Balena, and we will use anaconda and the intel compilers. So load up anaconda and intel:

	module load anaconda3/2.5.0
	module load intel/mpi/64/5.1.3.210

change into your `pimaim_mpi` directory and:

	make

You should now have a `pimaim_mpi` executable in this directory. Next clone ppfit, and switch to the development branch.

	git clone https://github.com/bjmorgan/ppfit
	cd ppfit
	git checkout develop
	
and add this directory to your pythonpath:

	export PYTHONPATH=/YOUR/PATH/TO/PPFIT/FOLDER/HERE/ppfit/
	
In the ppfit folder are two folders *ppfit/*, containing the fitting code, and *fitabinitio/* containing (click for more info):


<details><summary> <b>fitabinitio.py</b>: an example python script which sets up and runs the fitting </summary><p>

####  **`fitabinitio.py`**

```python
options_read = read_options( 'options.yml' )

run_options = Options(options_read)
```
Reads the options from `option.yml` and initialises an instance of the `Options` class run_options.

```python
comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
```

sets up the MPI communicator and gets the rank.

```python
if rank ==0:
 outfile = open('OUTPUT','w')

 config_read = read_options('configs.yml')
 configs=[ Configuration.from_dict( run_options,config) for config in sorted(config_read.values(),key=lambda x:sorted(x.keys()))]
else:
 configs = None


configs = comm.bcast(configs, root=0)
```
Opens the output file `OUTPUT` on the master node , reads the configurations from  `config.yml` and sets up an instance of the `Configuration` class for each configuration in the training set. Broadcasts these configs to all nodes.

```python
training_set = Training_Set(configs, run_options )

fitting_parameters = Fitting_Parameter_Set.from_parameters_file( 'PARAMS' )
potential_file = Potential_File( 'template_BaTiO3', fitting_parameters )

chi_squared_scaling = { 'forces':   options_read[ 'scaling' ][ 'forces' ],
                        'dipoles':  options_read[ 'scaling' ][ 'dipoles' ],
                        'stresses': options_read[ 'scaling' ][ 'stresses' ] }
```

Sets up an instance of  the `Training_set` class with the configs, gets an instance of  `Fitting_Parameter_Set` from the `PARAMS` file, initialises an instance of `Potential_file` by fleshing out the template file with the fitting paramters.

```python

chi_squared_scaling = { 'forces':   options_read[ 'scaling' ][ 'forces' ],
                        'dipoles':  options_read[ 'scaling' ][ 'dipoles' ],
                        'stresses': options_read[ 'scaling' ][ 'stresses' ] }
```

Simply sets the scaling of the chi_squared funtions for the forces, dipoles and stresses - giving each a relative weight scaling when included in the total sum of chi squared used by the optimiser. 

```python
optimise( sum_of_chi.evaluate, fitting_parameters, options_read )

MPI.COMM_WORLD.Barrier()
MPI.Finalize()
```

Finally does the optimisation, and finalises the MPI ranks. 
</p>
</details>

<details><summary>  <b>options.yml</b>: the options file in .yaml format</summary><p>

####  **`options.yml`**
The example `options.yml` file ([.yaml](http://yaml.org/about.html)  formatting) has a few sections (some of which are commented out), but they can be broadly spllit into two groups. The options which control the run from the `ppfit` side of things, and the options which control the optimisation.

First controlling the `ppfit` side of things:

In `mpi4py` mode ppfit will run a single master process, and several worker processes. At each step in the optimisation the configurations in the training set are added to a queue and sent to the worker processes until all configurations have been run. These configurations are run in parallel over the worker processes. Also each worker process runs the `pimaim_mpi` code in parallel. So there are two layers of parallelism. 

The first few lines tell ppfit the MD code being used is pimaim, where to find the executable and that the MD code will run in parallel using MPI:

```
calculation:
    code: pimaim
    exec: ~/bin/pimaim_mpi/pimaim_mpi
    code_mpi: True
```
The next line tells ppfit to run the configurations in parallel using mpi4py.

```
run_configs: mpi4py #serial, pool or mpi4py
```

The next few lines control the lower level of parallelism, i.e. the running of the `pimaim_mpi` executable, and tell ppfit to use the `mpiexec.hydra` executable to run each instance of `pimaim_mpi` on the configuration. `-np` is the string flag fed to mpiexec.hydra before the number of processes requested (this will be dependent on the MPI library you use, so if you don't use the same check). `exec_proc` tells each worker the number of processes to run `pimaim_mpi`on. `mpi_options` contains any options flags needed by the mpirun executable. Note the `-f` here: this is necessary and should come last - it is the machine name flag used by the mpiexec executable which tell mpiexec the name of the machines available and the number of cpus on each.  This will vary by library (Check the man page if not using the same MPI flavour).

```
mpi_exec: mpiexec.hydra
mpi_np: -np          # string for number of processes
exec_proc: 2    # execute each instance of code over exec_proc processes
 mpi_options: -genv I_MPI_FABRICS=dapl,ofa,tcp,tmi,ofi --rr --map-by --node --bind-to node -f
```

The next lines tell ppfit to run over 4 nodes, that there are 16 cpus on each node and that the total number of workers required is 32.

```
no_nodes:  4             # number of nodes
no_cpu_node: 16          # number of cpus per node
no_slaves: 32            # number of workers. Total no of mpi proc to request in submission is no_slaves + 1(master)
```

In the example above ppfit is run on 4 nodes, each containing 16 cpus, on each of these nodes 8 worker processes are started which each run `pimaim_mpi` on two processes. A couple of notes:

- The total number of workers should be less than the number of configurations, otherwise there will processes idle.
- When submitting the run script to the Balena queue the number of processes `fitabinitio.py` is run on should be `no_slaves+1`, as the master process needs to be included. If this number is incorrect the calculation will hang.
- The most efficient set-up will depend on the number of configurations, and the trainning set. It's recommended to do some testing.
</p>
</details>

<details><summary><b>run_fit_mpi4py_example.sh</b>: an example submission script to run on Balena using mpi4py</summary><p>

The `run_fit_mpi4py_example.sh` script is used to request the resources for your fitting calculation, and submit the job to the queue on Balena. There is more information on running jobs on Balena in the wiki [Balena wiki](https://wiki.bath.ac.uk/display/BalenaHPC), but there are a few non standard things in this run script, including various MPI flags that should be included (which will also vary by library used, so if you switch from the setup described you will need to check these). 

The most notable thing you might have to concern yourself with are the lines:

```
srun hostname -s | sort -u > slurm.hosts
echo $SLURM_NNODES

for ((i=0; i<$SLURM_NNODES; i++)); do
  echo ":16">>proc
done
paste slurm.hosts proc > temp
awk '{print $1$2}' temp > slurm.hosts

rm temp
rm proc
```

which generate a machine file called `slurm.host`, which tells MPI on which machine to run the calculations. If you run in interactive mode (using `sinteractive`), which can be convenient for testing the set-up of the calculation and to ensure that it is running smoothly, you will need to generate this file by hand. If this machine file is incorrect the calculation will hang.

Again note the number of processes you run `fitabinitio.py` on should be one more than the number of slaves requested in the options.yml file. 

</p>
</details>

<details><summary><b> run_fit_serial_pool_example.sh</b>: an example submission script to run in serial, or in multiprocessing mode on Balena</summary><p>
Simple submission script for serial and pool version. No need for a host file, or to export the mpi environment variables. The anaconda and intel modules should be loaded, and the `PYTHONPATH` environment variable exported though. 
</p>
</details>

<details><summary><b>configs.yml</b>: a file containing the configurations in the training set</summary><p>

For each member of your training set you should have an entry in the `configs.yml` file. The form of each entry is:

```
SEEDNAME:
    species:
        - Element 1
        - Element 2
        - Element 3
    directory:     configs/SEEDNAME/
    runtime_file:  runtime_SEEDNAME.inpt
    restart_file:  restart_SEEDNAME.dat
    forces_file:   SEEDNAME.force
    dipoles_file:  SEEDNAME.dip
    stresses_file: SEEDNAME.stress
```
where `SEEDNAME` is the string you have chosen to describe your member and the name of the folder in which to find the training data for this member. 

I've included a simple one-liner script [set_configs.sh](../blob/master/ppfit%20guide%20files/configs_setup/set_configs.sh) and a [base file](../blob/master/ppfit%20guide%20files/configs_setup/configs.base) for setting this `configs.yml` file up from a file including a list of your seednames line by line, `configs`. Edit the configs.base file to include the elements in your training sets and run this script as 

```
./set_configs.sh configs
```
to produce the yml file.


</p>
</details>

<details><summary><b> template_BaTiO3</b>: an example potential template</summary><p>
The template file is the skeleton of the ``potential.inpt`` file used during the MD run. This is the potential file we wish to optimise during the fitting. During each step of the fitting optimisation the parameters being optimised are substituted into this file and the errors between the ab-initio forces, dipoles and stresses and those produced during the MD are then optimised. This proceeds until the errors reach the required level and you have a fitted potential. 

First a bit about the potential file. This potential template, and probably the one you will most likely be interested in is a **DIP**ole **P**olarisable **I**on **M**odel (DIPPIM) potential ([1](https://www.sciencedirect.com/science/article/pii/S0166128006001412) ,[2](https://www.sciencedirect.com/science/article/pii/S0031920107000623),[3](http://dx.doi.org/10.1088/0953-8984/23/25/255402)). There is no documentation for the pimaim code, and I have found the best way to figure out the format for the potentials is to look through the source. Th e relevant part of the code to look at is ``readin.f90``, and the reading of the potential paramters starts around line 230.

Briefly we will take a look at the DIPPIM potential format. DIPPIM is made up of three components: short range repulsion, dispersion and polarisability. 

![](https://github.com/bjmorgan/group_info/blob/master/ppfit%20guide%20files/Equations/CodeCogsEqn1.png)

The short range repulsive interactions are given by:

![](https://github.com/bjmorgan/group_info/blob/master/ppfit%20guide%20files/Equations/CodeCogsEqn2.png)

The dispersion interactions account for the dipole-dipole and dipole-quadrupole terms:

![](https://github.com/bjmorgan/group_info/blob/master/ppfit%20guide%20files/Equations/CodeCogsEqn3.png)

The polarization interactions are represented by:

![](https://github.com/bjmorgan/group_info/blob/master/ppfit%20guide%20files/Equations/CodeCogsEqn4.png)

where q represents the charge, &theta; represents the x,y or z cartesian axes, &mu;<sub>&theta;</sub><sup>*i*</sup>  is the &theta; component of the dipole moment if ion *i*, &Delta;<sub>*ij*</sub>=1 if *i* and *j* refer to different ionic species and zero otherwise, &alpha;<sup>*i*</sup>  is the dipole polarizability of ion *i* and *k*<sup>*i*</sup><sub>1</sub> = 1/(2&alpha;<sup>*i*</sup>). T<sub>&theta;</sub><sup>(1)</sup> and *T*<sup>(2)</sup><sub>&theta;<sub>1</sub>&theta;<sub>2</sub></sub> are the charge-dipole and dipole-dipole stress tensors. 

*f*(*r*<sub>*ij*</sub>) are the Tang-Toennies dispersion damping function [4](https://aip.scitation.org/doi/abs/10.1063/1.447150):

![](https://github.com/bjmorgan/group_info/blob/master/ppfit%20guide%20files/Equations/CodeCogsEqn5.png)

where the *b*<sup>*ij*</sup> term determines the range of the damping and the *c*<sup>*ij*</sup> term determines the strength of the response.


Below is the potential format for a XFT DIPPIM potential Li-Al doped magnesium spinel (Li<sub>*x*</sub>Mg<sub>1-2*x*</sub>Al<sub>2+*x*</sub>O<sub>4</sub>). The atom ordering is O-Mg-Li-Al. 

![](https://github.com/bjmorgan/group_info/blob/master/ppfit%20guide%20files/Images/potential.png)
</p>
</details>

<details><summary><b> PARAMS</b>: the starting parameters for your potential to flesh out the template file</summary><p>

</p>
</details>






	
	















