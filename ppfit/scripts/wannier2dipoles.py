#! /usr/bin/env python3
'''
This script calculates dipoles associated with atoms from wannier cantres, and outputs a file
containing forces, a file containing stresses and a restart file all for use with ppfit

Usage:
 'wannier2dipoles.py' 
Command line flags:
  -h, --help: help statement and exit
Command line arguments (Optional):
  -p,--poscar  : Structure file in VASP POSCAR format. Default = 'POSCAR'
  -wc, -wannier: File containing wannier centres in xyz format. Default = 'wannier90_centres.xyz'
  -u,--units   : "Output units: either 'Bohr' or 'Angs'. Input is assumed Angs. Defualt = 'Angs'
  -rc, --rcut  : Cut-off radius for wannier centres to be included in dipole calc. Default = 0.7 Angs
  -s, --sort   : Reorder the atoms? Default is false.
  -o, -order   : If reordering what is the new atomic order (must match atomic species in POSCAR).
  -d, -dipoles : get the dipoles if there is a wannier90_centres.xyz. Default is False
Requires:
  pymatgen, elementtree and numpy

'''
import sys
import numpy as np
import math
import argparse 
from pymatgen import Lattice, Structure, Molecule
from pymatgen.io.cif import CifParser
import xml.etree.ElementTree as ET

def parse_command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--poscar",default='POSCAR',help="Structure file in VASP POSCAR format. Default = 'POSCAR'")
    parser.add_argument("-wf","--wannier",default='wannier90_centres.xyz',help="File containing wannier centres in xyz format. Default = 'wannier90_centres.xyz'")
    parser.add_argument("-u","--units",choices = ['Bohr', 'Angs'], default = 'Angs', help="Output units: either 'Bohr' or 'Angs'. Input is assumed Angs. Defualt = 'Angs'")
    parser.add_argument("-rc","--rcut",type=float,  default = 0.7, help="Cut-off radius for wannier centres to be included in dipole calc. Default = 0.7 Angs" )
    parser.add_argument("-s","--sort",choices = ['True', 'False'], default=False, help="Would you like to resort the order of the elements when outputting?. Default = False")
    parser.add_argument("-o","--order",nargs='+',default=None, help="If -s (--sort) == True -o will give the species order to resort the output. Default = None")
    parser.add_argument("-d","--dipoles",nargs='+',default=False, help="If -d == True the wannier_centres.xyz file will be used to calculate the dipoles, otherwise just the forces and stresses are output. Default = False")
    args = parser.parse_args()
    return args

def main():
    args = parse_command_line_arguments()
    if args.units == 'Angs':
        fac = -2
    elif args.units=='Bohr':
        fac= -2*1.889725989
    else:
       print('Error: units should be either \'Bohr\' or \'Angs\'. \'pymatgen_dipoles.py -h\' for more details')
       exit()
    #Parse vasprun for the forces and stresses
    tree = ET.parse('vasprun.xml')
    root = tree.getroot()
    stresses = np.array(
               [ v.text.split() for v in root.findall('.//calculation/varray[@name="stress"]/v')]
               ,dtype='float')

    forces=[]
    for v in root.findall('.//calculation/varray[@name="forces"]/v'):
        forces.append([float(v.text.split()[0]),float(v.text.split()[1]),float(v.text.split()[2])])

    # set up pymatgen structure object from POSCAR and change ordering or it and forces if instructed

    structure_pos = Structure.from_file(args.poscar)

    if (args.sort):
       symbols=[species for species in structure_pos.symbol_set]
       if (set(symbols) == set(args.order)):
           structure_sorted=Structure(lattice=structure_pos.lattice,species=[],coords=[])
           structure_restart=Structure(lattice=structure_pos.lattice,species=[],coords=[])
           forces_sorted=[]
           for symbol in args.order:
              for i,site in enumerate(structure_pos.sites):
                  if (site.species_string==symbol):
                     structure_sorted.append(symbol,site.coords,coords_are_cartesian=True)
                     structure_restart.append(symbol,site.coords,coords_are_cartesian=True)
                     forces_sorted.append(forces[i])
           structure_pos=structure_sorted
           forces=forces_sorted
       else:
           print('Error: elements in list passed with "-o(--order)" does not match that found in POSCAR')
           print('Passed: ',args.order)
           print('POSCAR: ',symbols)
           exit()
    else:
        structure_restart = Structure.from_file(args.poscar)
    
    
    natoms = len(structure_pos.sites)
    ## set up pymatgen molecule for wannier centres and transform to POSCAR lattice, then append to structure
    if (args.dipoles):
        wc = Molecule.from_file(args.wannier)
    
        # Create new structure containing the WC and atoms in the POSCAR lattice 
        coords = np.dot(wc.cart_coords,structure_pos.lattice.inv_matrix)
        symbols = wc.species
        ncentres = len(wc.sites)
    
        [structure_pos.append(symbol,coord) for symbol,coord in zip(symbols,coords)]
    
        structure_pos.to(filename="atomsNcentres.cif")
    
        dip_file = open('dipoles.out','w')
        centres_file = open('centres.out','w')
    
        # Calculate the dipoles:
        index = 0
        wc_attributed = 0
        for atom in range(natoms):
            index += 1
            dx=0
            dy=0
            dz=0  
            nwann = 0
            for wan in range(natoms,natoms+ncentres): 
                r = structure_pos.sites[atom].distance_and_image(structure_pos.sites[wan])
    
                if r[0] < args.rcut:
                    nwann += 1
                    image_pos = structure_pos.sites[wan].frac_coords+r[1]
                    dx = dx+structure_pos.lattice.get_cartesian_coords(image_pos)[0] - structure_pos.sites[atom].coords[0] 
                    dy = dy+structure_pos.lattice.get_cartesian_coords(image_pos)[1] - structure_pos.sites[atom].coords[1] 
                    dz = dz+structure_pos.lattice.get_cartesian_coords(image_pos)[2] - structure_pos.sites[atom].coords[2]
            wc_attributed = wc_attributed + nwann
            centres_file.write('{0:5d} {1:5d} \n'.format(index,nwann))
            dip_file.write('{0:5d} {1:< .7f} {2:< .7f} {3:< .7f}\n'.format(index,fac*dx,fac*dy,fac*dz))
    
        if wc_attributed != ncentres:
            print('I found '+str(wc_attributed)+' centres out of '+str(ncentres)+'\n')
            print('Please change the cutoff or check your input structures\n')
    
    
    #  write out restart file
    
    restart_file= open('restart.dat','w')
    
    restart_file.write('T\n')
    restart_file.write('F\n')
    restart_file.write('F\n')
    restart_file.write('F\n')
    
    if args.units == 'Angs':
      for site in structure_restart.sites:
          restart_file.write('{0:< .7f} {1:< .7f} {2:< .7f}\n'.format(site.frac_coords[0]*structure_restart.lattice.a,site.frac_coords[1]*structure_restart.lattice.b,site.frac_coords[2]*structure_restart.lattice.c))
    else:
      for site in structure_restart.sites:
          restart_file.write('{0:< .7f} {1:< .7f} {2:< .7f}\n'.format(structure_restart.lattice.a*site.frac_coords[0]*fac/-2,structure_restart.lattice.b*site.frac_coords[1]*fac/-2,structure_restart.lattice.c*site.frac_coords[2]*fac/-2))
    
    
    count=0
    for line in structure_restart.lattice.matrix.transpose():
       restart_file.write('{0:< .7f} {1:< .7f} {2:< .7f}\n'.format(line[0]/structure_restart.lattice._lengths[count],line[1]/structure_restart.lattice._lengths[count],line[2]/structure_restart.lattice._lengths[count]))
    
       count+=1
    
    count=0
    if args.units == 'Angs':
        restart_file.write('  {0:^ .9f} \n'.format(structure_restart.lattice.a))
        restart_file.write('  {0:^ .9f} \n'.format(structure_restart.lattice.b))
        restart_file.write('  {0:^ .9f} \n'.format(structure_restart.lattice.c))
    else:
        restart_file.write('  {0:^ .9f} \n'.format(structure_restart.lattice.a*1.889725989))
        restart_file.write('  {0:^ .9f} \n'.format(structure_restart.lattice.b*1.889725989))
        restart_file.write('  {0:^ .9f} \n'.format(structure_restart.lattice.c*1.889725989))
    
    
    
    # write out stress file
    stresses_file=open('stresses.out','w')
    
    if args.units == 'Angs':
       stresses_file.write('{0:< .7f} {1:< .7f} {2:< .7f}\n'.format(stresses[0,0],stresses[1,1],stresses[2,2]))
       stresses_file.write('{0:< .7f} {1:< .7f} {2:< .7f}\n'.format(stresses[0,1],stresses[0,2],stresses[1,2]))
    else:
       stresses_file.write('{0:< .7f} {1:< .7f} {2:< .7f}\n'.format(stresses[0,0]/294219.1219,stresses[1,1]/294219.1219,stresses[2,2]/294219.1219))
       stresses_file.write('{0:< .7f} {1:< .7f} {2:< .7f}\n'.format(stresses[0,1]/294219.1219,stresses[0,2]/294219.1219,stresses[1,2]/294219.1219))
    
    #write out forces file
    forces_file=open('forces.out','w')
    
    if args.units == 'Angs':
       for f in forces:
           forces_file.write('{0:^ .9f} {1:^ .9f} {2:^ .9f}\n'.format(f[0],f[1],f[2]))
    else:
       for f in forces:
           forces_file.write('{0:^ .9f} {1:^ .9f} {2:^ .9f}\n'.format(f[0]/51.42208619083232,f[1]/51.42208619083232,f[2]/51.42208619083232))

if __name__ == '__main__':
    main()
