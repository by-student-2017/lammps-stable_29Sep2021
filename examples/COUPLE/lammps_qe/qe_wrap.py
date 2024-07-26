#!/usr/bin/env python

# ----------------------------------------------------------------------
# LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
# http://lammps.sandia.gov, Sandia National Laboratories
# Steve Plimpton, sjplimp@sandia.gov
# ----------------------------------------------------------------------

# Syntax: vasp_wrap.py file/zmq POSCARfile

# wrapper on VASP to act as server program using CSlib
#   receives message with list of coords from client
#   creates VASP inputs
#   invokes VASP to calculate self-consistent energy of that config
#   reads VASP outputs
#   sends message with energy, forces, pressure to client

# NOTES:
# check to insure basic VASP input files are in place?
# could archive VASP input/output in special filenames or dirs?
# need to check that POTCAR file is consistent with atom ordering?
# could make syntax for launching VASP more flexible
#   e.g. command-line arg for # of procs
# detect if VASP had an error and return ERROR field, e.g. non-convergence ??

from __future__ import print_function
import sys
import re

version = sys.version_info[0]
if version == 3:
  sys.exit("The CSlib python wrapper does not yet support python 3")

import subprocess
import xml.etree.ElementTree as ET
from cslib import CSlib

# comment out 2nd line once 1st line is correct for your system

#vaspcmd = "srun -N 1 --ntasks-per-node=1 " + \
#          "-n 1 dftb+"
#vaspcmd = "touch tmp"
vaspcmd = "mpirun -np 8 /usr/bin/pw.x < scf.in > scf.out"

# enums matching FixClientMD class in LAMMPS

SETUP,STEP = range(1,2+1)
DIM,PERIODICITY,ORIGIN,BOX,NATOMS,NTYPES,TYPES,COORDS,UNITS,CHARGE = range(1,10+1)
FORCES,ENERGY,VIRIAL,CHARGES,ERROR = range(1,5+1)

# -------------------------------------
# functions

# error message and exit

def error(txt):
  print("ERROR:",txt)
  sys.exit(1)

# -------------------------------------
# read initial VASP POSCAR file to setup problem
# return natoms,ntypes,box

def vasp_setup(poscar): # poscar = gen format

  qe_input = open(poscar,'r')
  
  box = []
  types_symbol = []
  QMatoms = 0
  ntypes = 0
  
  while True:
    line = qe_input.readline()
    words = line.split()
    if not line: break
    
    #print(words)
    
    if words[0] == 'conv_thr':
      print("\n Note: setting 1 meV/atom \n")
      print("   conv_thr = "+str(1.0e-3/13.602*QMatoms)+" [Ry] \n")
    
    # pattern match for nat (QMatoms)
    if words[0] == 'nat':
      QMatoms = int(words[2])      #print(QMatoms)
    elif words[0] == 'nat=':
      QMatoms = int(words[1])      #print(QMatoms)
    
    # pattern match for ntyp (ntypes)
    if words[0] == 'ntyp':
      ntypes = int(words[2])
    elif words[0] == 'ntyp=':
      ntypes = int(words[1])
    
    # box size
    # The lattice parameter alat = A (in ANGSTROM ).
    if words[0] == 'A':
      A = float(words[2])
    elif words[0] == 'A=':
      A = float(words[2])
    #
    if words[0] == 'CELL_PARAMETERS':
      for i in range(3):
        line = qe_input.readline()
        vectors = line.strip().split()
        box += [float(value)*A for value in vectors]
    
    # types_symbol
    if words[0] == 'ATOMIC_POSITIONS':
      for i in range(QMatoms):
        line = qe_input.readline()
        words = line.strip().split()
        types_symbol += [words[0]]
  
  #print("nat:",QMatoms)
  #print("ntyp:",ntypes)
  #print("The lattice parameter alat = A (in ANGSTROM ):",A)
  #print("lattic vectors: ",box)
  #print("atomic symbol of types:",types_symbol)
  
  return QMatoms,ntypes,box,A,types_symbol
  
# -------------------------------------
# write a new POSCAR file for VASP

def poscar_write(poscar,QMatoms,ntypes,types,coords,box,A,types_symbol):

  template = open(poscar,'r')
  new_input = open("scf.in",'w')
  
  box_block = False; i = 0
  geometry_block = False; j = 0
  
  while True:
    line = template.readline()
    if not line: break
    
    if re.search("conv_thr",line):
      line = "conv_thr = "+str(1.0e-3/13.602*QMatoms)+" , \n" # 1 meV/atom
    
    if box_block and i < 3:
      x = box[3*i+0]/A
      y = box[3*i+1]/A
      z = box[3*i+2]/A
      coord_string = " %g %g %g \n" % (x,y,z)
      line = coord_string
      i += 1
    
    if re.search("CELL_PARAMETERS",line):
      box_block = True
    
    if geometry_block and j < QMatoms:
      t = types_symbol[j]
      x = coords[3*j+0]/A
      y = coords[3*j+1]/A
      z = coords[3*j+2]/A
      coord_string = "%s  %g %g %g \n" % (t,x,y,z)
      line = coord_string
      j += 1
    
    if re.search("ATOMIC_POSITIONS",line):
      geometry_block = True
    
    print(line,file=new_input,end='')

  new_input.close()

# -------------------------------------
# read a VASP output vasprun.xml file
# uses ElementTree module
# see https://docs.python.org/2/library/xml.etree.elementtree.html

def vasprun_read():
  tree = ET.parse('pwscf.xml')
  root = tree.getroot()
  
  fout = []
  sout = []
  qout = []
  
  #---------------------------------------------------------------------------------------
  # https://courses.grainger.illinois.edu/mse404ela/sp2019/6.DFT-walkthrough.html
  # Quantum Espresso uses "atomic units" (Ryd for energy, Bohr for distance). 
  # You can convert these to units you will find more convenient:
  # 1 bohr = 1 a.u. (atomic unit) = 0.529177211 Angstrom
  # 1 Rydberg (Ryd or Ry) = 13.60569301 eV = 0.5 Ha
  # 1 kcal/mol = 43.36 meV
  # 1 kJ/mol = 10.36 meV
  # 1 eV = 1.602176621e-19 J
  # 1 eV/Angstrom^3 = 160.2176621 GPa
  # 1 Ryd/bohr^3 = 14,710.5076 GPa
  # 1 Ha/bohr^3 = 2 Ry/bohr^3 = 2*14710.5076 GPa = a.u. = au
  # 1 bar = 100 kPa
  #---------------------------------------------------------------------------------------
  # pwscf.xml <!--All quantities are in Hartree atomic units unless otherwise specified-->
  #---------------------------------------------------------------------------------------
  # Not "Rydberg (Ry) atomic units (screen)", but "Hartree (Ha) atomic units (xml file)"
  #---------------------------------------------------------------------------------------
  for output in root.findall('output'):
    for neighbor in output.iter('atomic_structure'):
      QMatoms = int(neighbor.attrib.get('nat'))
    
    eunitconv=1.0/0.367493245336341E-01 # Ha => eV
    for total_energy in output.iter('total_energy'):
      eout = float(total_energy.find('etot').text)*eunitconv # eV
    
    funitconv=1.0/0.194469064593167E-01 # Ha/Bohr = Ha/au => eV/Angstrom
    forces = output.find('forces').text
    fxyz = forces.split()
    fout = [float(value)*funitconv for value in fxyz]
    
    sunitconv=1.0/0.339893208050290E-13 # Ha/bohr^3 = 2 Ry/bohr^3 = au => Pa
    stress = output.find('stress').text
    sxyz = stress.split()
    stensor = [float(value)*sunitconv for value in sxyz]
  
  # 0 1 2 = [0][0] [0][1] [0][2]
  # 3 4 5 = [1][0] [1][1] [1][2]
  # 6 7 8 = [2][0] [2][1] [2][2]
  sxx = stensor[0] # sxx = stensor[0][0]
  syy = stensor[4] # syy = stensor[1][1]
  szz = stensor[8] # szz = stensor[2][2]
  # symmetrize off-diagonal components
  sxy = 0.5 * (stensor[1] + stensor[3]) # sxy = 0.5 * (stensor[0][1] + stensor[1][0])
  sxz = 0.5 * (stensor[2] + stensor[6]) # sxz = 0.5 * (stensor[0][2] + stensor[2][0])
  syz = 0.5 * (stensor[5] + stensor[7]) # syz = 0.5 * (stensor[1][2] + stensor[2][1])
  sout = [sxx,syy,szz,sxy,sxz,syz]
  
  # dummy
  #for i in range(QMatoms):
  #  qout += [float(0.0)]
    
  #print(QMatoms)
  #print(eout)
  #print(fout)
  #print(sout)
  
  return eout,fout,sout,qout

# -------------------------------------
# main program

# command-line args

if len(sys.argv) != 3:
  print("Syntax: python vasp_wrap.py file/zmq POSCARfile")
  sys.exit(1)

mode = sys.argv[1]
poscar_template = sys.argv[2]

if mode == "file": cs = CSlib(1,mode,"tmp.couple",None)
elif mode == "zmq": cs = CSlib(1,mode,"*:8888",None)
else:
  print("Syntax: python vasp_wrap.py file/zmq POSCARfile")
  sys.exit(1)

QMatoms = 0
types_symbol = []
QMatoms,ntypes,box,A,types_symbol = vasp_setup(poscar_template)

# initial message for MD protocol

msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
if msgID != 0: error("Bad initial client/server handshake")
protocol = cs.unpack_string(1)
if protocol != "md": error("Mismatch in client/server protocol")
cs.send(0,0)

# endless server loop

while 1:

  # recv message from client
  # msgID = 0 = all-done message

  msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
  if msgID < 0: break

  # SETUP receive at beginning of each run
  # required fields: DIM, PERIODICTY, ORIGIN, BOX, 
  #                  NATOMS, NTYPES, TYPES, COORDS
  # optional fields: others in enum above, but VASP ignores them

  if msgID == SETUP:
    
    origin = []
    box = []
    natoms_recv = ntypes_recv = 0
    types = []
    coords = []
    lmp_charges = []
    
    for field in fieldID:
      if field == DIM:
        dim = cs.unpack_int(DIM)
        if dim != 3: error("VASP only performs 3d simulations")
      elif field == PERIODICITY:
        periodicity = cs.unpack(PERIODICITY,1)
        #if not periodicity[0] or not periodicity[1] or not periodicity[2]:
        #  error("VASP wrapper only currently supports fully periodic systems")
      elif field == ORIGIN:
        origin = cs.unpack(ORIGIN,1)
      elif field == BOX:
        box = cs.unpack(BOX,1)
      elif field == NATOMS:
        natoms_recv = cs.unpack_int(NATOMS)
        natoms = natoms_recv
        #if natoms != natoms_recv:
        #  error("VASP wrapper mis-match in number of atoms")
      elif field == NTYPES:
        ntypes_recv = cs.unpack_int(NTYPES)
        #if ntypes != ntypes_recv:
        #  error("VASP wrapper mis-match in number of atom types")
      elif field == TYPES:
        types = cs.unpack(TYPES,1)
      elif field == COORDS:
        coords = cs.unpack(COORDS,1)
      elif field == CHARGE:
        lmp_charges = cs.unpack(CHARGE,1)

    if not origin or not box or not natoms or not ntypes or \
       not types or not coords:
      error("Required VASP wrapper setup field not received");

  # STEP receive at each timestep of run or minimization
  # required fields: COORDS
  # optional fields: ORIGIN, BOX

  elif msgID == STEP:

    coords = []
    
    for field in fieldID:
      if field == COORDS:
        coords = cs.unpack(COORDS,1)
      elif field == ORIGIN:
        origin = cs.unpack(ORIGIN,1)
      elif field == BOX:
        box = cs.unpack(BOX,1)
    
    if not coords: error("Required VASP wrapper step field not received");

  else: error("VASP wrapper received unrecognized message")
  
  # create POSCAR file
  
  poscar_write(poscar_template,QMatoms,ntypes,types,coords,box,A,types_symbol)

  # invoke DFTB+
  
  #print("\nLaunching DFTB+ ...")
  #print(vaspcmd)
  subprocess.check_output(vaspcmd,stderr=subprocess.STDOUT,shell=True)
  
  # process VASP output
    
  energy,forces,virial,charges = vasprun_read()
  
  #for i in range(QMatoms,natoms):
  #  charges += [float(lmp_charges[i])]
  
  for i in range(natoms):
    charges += [float(lmp_charges[i])]
  
  # debag
  #print(forces)
  #print(virial)
  #print(charges)

  # convert VASP kilobars to bars

  for i,value in enumerate(virial): virial[i] *= 1.0/(100.0*1000.0) # Pa to bars (1 bar = 0.1 MPa = 100 kPa)

  # return forces, energy, pressure, charge to client
  
  cs.send(msgID,4);
  cs.pack(FORCES,4,3*natoms,forces)
  cs.pack_double(ENERGY,energy)
  cs.pack(VIRIAL,4,6,virial)
  cs.pack(CHARGES,4,natoms,charges)
  
# final reply to client
  
cs.send(0,0)

# clean-up

del cs
