#!/usr/bin/env python

# ----------------------------------------------------------------------
# LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
# http://lammps.sandia.gov, Sandia National Laboratories
# Steve Plimpton, sjplimp@sandia.gov
# ----------------------------------------------------------------------

# Syntax: nwchem_wrap.py file/zmq ao/pw input_template
# file/zmg = messaging mode, must match LAMMPS messaging mode
# ao/pw = basis set mode, selects between atom-centered and plane-wave
#         the input_template file must correspond to the appropriate basis set mode:
#         the "ao" mode supports the scf and dft modules in NWChem,
#         the "pw" mode supports the nwpw module.
# input_template = NWChem input file used as template, must include a
#                  "geometry" block with the atoms in the simulation, dummy
#                  xyz coordinates should be included (but are not used).
#                  Atom ordering must match LAMMPS input.

# wrapper on NWChem
#   receives message with list of coords
#   creates NWChem inputs
#   invokes NWChem to calculate self-consistent energy of that config
#   reads NWChem outputs
#   sends message with energy, forces, pressure to client

from __future__ import print_function
import sys

version = sys.version_info[0]
if version == 3:
  sys.exit("The CSlib python wrapper does not yet support python 3")

import subprocess
import re
import os
import shutil # for copy command
from cslib import CSlib

# comment out 2nd line once 1st line is correct for your system

cp2kcmd = "mpirun -np 1 /usr/bin/cp2k -i cp2k.inp -o cp2k.out"
#cp2kcmd = "touch tmp"

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
# read initial input file to setup problem
# return natoms

def nwchem_setup_ao(input):

  template = open(input,'r')

  geometry_block = False
  natoms = 0

  line = template.readline()
  natoms = int(line.split()[0])

  #print(natoms)
  
  return natoms

# -------------------------------------
# write a new input file for NWChem
# assumes the NWChem input geometry is to be specified in angstroms

def nwchem_input_write_ao(input,coords):

  template = open(input,'r')
  new_input = open("cp2k.xyz",'w')
  ixyz = 0
  j = 0

  while True:
    line = template.readline()
    if not line: break

    ixyz += 1
    if ixyz == 1:
      pass
      natoms = line.split()[0]
    if ixyz == 2:
      pass
    if ixyz >= 3:
      x = coords[3*j+0]
      y = coords[3*j+1]
      z = coords[3*j+2]
      coord_string = "  %g %g %g \n" % (x,y,z)
      atom_string = line.split()[0]
      line = atom_string + coord_string
      j += 1

    print(line,file=new_input,end='')

  new_input.close()

# -------------------------------------
# read a NWChem output nwchem_lammps.out file

def nwchem_read_ao(natoms, log):
  cp2k_output = open(log, 'r')

  eout = 0.0
  fout = []
  reading_forces = False
  qout = []
  reading_charges = False
  reading_charges_flag = False
  forces_pattern = r"^\s*\d+\s+\d+\s*\w{1,2}(?:\s+-?(?:\d+.?\d*|\d*.?\d+)){3}"

  lines = cp2k_output.readlines()
  
  start_index = None
  for i, line in enumerate(reversed(lines)):
    #if ' TOTAL NUMBERS AND MAXIMUM NUMBERS' in line:
    if re.search("TOTAL NUMBERS AND MAXIMUM NUMBERS", line):
      start_index = len(lines) - i - 1
      #print(start_index)
      break

  if start_index is not None:
    for line in lines[start_index:]:
      if not line: break
      #print(line)

      # pattern match for energy
      #eunitconv=1.0/0.367493245336341E-01
      if re.search("Total energy:",line):
        eout = float(line.split()[2])

      # pattern match for forces
      #funitconv=1.0/0.194469064593167E-01
      if re.search("SUM OF ATOMIC FORCES", line):
        reading_forces = False
      if reading_forces and not re.match("#",line) and re.search(forces_pattern,line):
        forces = line.split()
        fout += [float(forces[3]), float(forces[4]), float(forces[5])]
      if re.search("ATOMIC FORCES in \[a.u.\]",line):
        reading_forces = True

      # pattern match for charges
      if re.search("# Total charge", line):
        reading_charges = False
        reading_charges_flag = False
      if reading_charges and reading_charges_flag:
        charges = line.split()
        #qout += [float(charges[2]), float(charges[3]), float(charges[4])]
        qout += [float(charges[4])]
      if re.search("Mulliken Population Analysis",line):
        reading_charges = True
      if re.search("Net charge",line):
        reading_charges_flag = True

  # convert units
  hartree2eV    = 27.21138602   # a.u. = Ha = 2 Ry = 2*13.602 eV
  bohr2angstrom = 0.52917721092 # 
  eout = eout * hartree2eV
  fout = [i * hartree2eV/bohr2angstrom for i in fout]

  #print(eout)
  #print(fout)
  #print(qout)

  return eout,fout,qout

# -------------------------------------
# read initial planewave input file to setup problem
# return natoms,box
def nwchem_setup_pw(input):

  template = open(input,'r')

  COORD_block = False
  coord_pattern = r"^\s*\w{1,2}(?:\s+-?(?:\d+.?\d*|\d*.?\d+)){3}"
  natoms = 0
  box = []

  while True:
    line = template.readline()
    if not line: break

    if re.search("\&CELL",line) and not re.match("#",line):
      for i in range(3):
        line = template.readline()
        box += [float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]


    if COORD_block and re.search("\&END COORD",line) and not re.match("#",line):
      COORD_block = False

    if COORD_block and not re.match("#",line) and re.search(coord_pattern,line):
      natoms += 1

    if not COORD_block and re.search("\&COORD",line) and not re.match("#",line):
      COORD_block = True

  #print(natoms)
  #print(box)

  return natoms,box

# -------------------------------------
# write a new planewave input file for NWChem
# assumes the NWChem input geometry is to be specified fractional coordinates

def nwchem_input_write_pw(input,coords,box):

  template = open(input,'r')
  new_input = open("cp2k.inp",'w')
  
  CELL_block = False
  COORD_block = False
  coord_pattern = r"^\s*\w{1,2}(?:\s+-?(?:\d+.?\d*|\d*.?\d+)){3}"
  i = 0 # for CELL
  j = 0 # for COORD

  while True:
    line = template.readline()
    if not line: break

    if CELL_block and re.search("\&END CELL",line) and not re.match("#",line):
      CELL_block = False

    if CELL_block and not re.match("#",line):
      x = box[3*i+0]
      y = box[3*i+1]
      z = box[3*i+2]
      cell_string = "  %g %g %g \n" % (x,y,z)
      atom_string = line.split()[0]
      line = "      "+atom_string + cell_string
      i += 1

    if not CELL_block and re.search("\&CELL",line) and not re.match("#",line):
      CELL_block = True


    if COORD_block and re.search("\&END COORD",line) and not re.match("#",line):
      COORD_block = False

    if COORD_block and not re.match("#",line) and re.search(coord_pattern,line):
      x = coords[3*j+0]
      y = coords[3*j+1]
      z = coords[3*j+2]
      coord_string = "  %g %g %g \n" % (x,y,z)
      atom_string = line.split()[0]
      line = "      "+atom_string + coord_string
      j += 1

    if not COORD_block and re.search("\&COORD",line) and not re.match("#",line):
      COORD_block = True

    print(line,file=new_input,end='')

  new_input.close()

# -------------------------------------
# read a NWChem output nwchem_lammps.out file for planewave calculation

def nwchem_read_pw(log):
  cp2k_output = open(log, 'r')

  eout = 0.0
  stensor = []
  sout = []
  reading_stress = False
  fout = []
  reading_forces = False
  qout = []
  reading_charges = False
  reading_charges_flag = False
  forces_pattern = r"^\s*\d+\s+\d+\s*\w{1,2}(?:\s+-?(?:\d+.?\d*|\d*.?\d+)){3}"
  istress = 0

  lines = cp2k_output.readlines()
  
  start_index = None
  for i, line in enumerate(reversed(lines)):
    #if ' TOTAL NUMBERS AND MAXIMUM NUMBERS' in line:
    if re.search("TOTAL NUMBERS AND MAXIMUM NUMBERS", line):
      start_index = len(lines) - i - 1
      #print(start_index)
      break

  if start_index is not None:
    for line in lines[start_index:]:
      if not line: break
      #print(line)

      # pattern match for energy
      #eunitconv=1.0/0.367493245336341E-01
      if re.search("Total energy:",line):
        eout = float(line.split()[2])

      # pattern match for forces
      #funitconv=1.0/0.194469064593167E-01
      if re.search("SUM OF ATOMIC FORCES", line):
        reading_forces = False
      if reading_forces and not re.match("#",line) and re.search(forces_pattern,line):
        forces = line.split()
        fout += [float(forces[3]), float(forces[4]), float(forces[5])]
      if re.search("ATOMIC FORCES in \[a.u.\]",line):
        reading_forces = True

      # pattern match for stress
      #sunitconv=1.0/0.367493245336341E-01
      if re.search("STRESS\|1\/3 Trace",line):
        reading_stress = False
      if reading_stress and istress < (3+1):
        if not istress == 0:
          stress = line.split()
          #print(stress)
          stensor += [float(stress[2]), float(stress[3]), float(stress[4])]
        istress += 1
      if re.search("Analytical stress tensor \[GPa\]",line):
        reading_stress = True

      # pattern match for charges
      if re.search("# Total charge", line):
        reading_charges = False
        reading_charges_flag = False
      if reading_charges and reading_charges_flag:
        charges = line.split()
        #qout += [float(charges[2]), float(charges[3]), float(charges[4])]
        qout += [float(charges[4])]
      if re.search("Mulliken Population Analysis",line):
        reading_charges = True
      if re.search("Net charge",line):
        reading_charges_flag = True

  sxx = stensor[0]
  syy = stensor[4]
  szz = stensor[8]
  sxy = 0.5 * (float(stensor[1]) + float(stensor[3]))
  sxz = 0.5 * (float(stensor[2]) + float(stensor[6]))
  syz = 0.5 * (float(stensor[5]) + float(stensor[7]))
  sout = [sxx,syy,szz,sxy,sxz,syz]

  # convert units
  hartree2eV    = 27.21138602   # a.u. = Ha = 2 Ry = 2*13.602 eV
  bohr2angstrom = 0.52917721092 # 
  #austress2bar  = 294210156.97  # bar = 0.1 kPa = 0.1 MPa = 0.1/1000 GPa
  gpa2bar       = 1.0e4
  eout = eout * hartree2eV
  fout = [i * hartree2eV/bohr2angstrom for i in fout]
  #sout = [i * austress2bar for i in sout]
  sout = [i * gpa2bar for i in sout]

  #print(eout)
  #print(fout)
  #print(sout)
  #print(qout)

  return eout,fout,sout,qout

# -------------------------------------
# main program

# command-line args
#
if len(sys.argv) != 4:
  print("Syntax: python nwchem_wrap.py file/zmq ao/pw input_template")
  sys.exit(1)

comm_mode = sys.argv[1]
basis_type = sys.argv[2]
input_template = sys.argv[3]

if comm_mode == "file": cs = CSlib(1,comm_mode,"tmp.couple",None)
elif comm_mode == "zmq": cs = CSlib(1,comm_mode,"*:5555",None)
else:
  print("Syntax: python nwchem_wrap.py file/zmq")
  sys.exit(1)


natoms = 0
box = []
if basis_type == "ao":
  natoms = nwchem_setup_ao(input_template)
elif basis_type == "pw":
  natoms,box = nwchem_setup_pw(input_template)

# initial message for AIMD protocol

msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
if msgID != 0: error("Bad initial client/server handshake")
protocol = cs.unpack_string(1)
if protocol != "md": error("Mismatch in client/server protocol")
cs.send(0,0)

# endless server loop

i = 0
directory_path = "cp2k_logs"
if not os.path.exists(directory_path):
  os.mkdir(directory_path)

file_path = 'cp2k.out'
if os.path.exists(file_path):
    os.remove(file_path)

while 1:

  # recv message from client
  # msgID = 0 = all-done message

  msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
  if msgID < 0: break

  # SETUP receive at beginning of each run
  # required fields: DIM, PERIODICITY, ORIGIN, BOX,
  #                  NATOMS, COORDS
  # optional fields: others in enum above, but NWChem ignores them

  if msgID == SETUP:

    origin = []
    box_lmp = []
    natoms_recv = ntypes_recv = 0
    types = []
    coords = []

    for field in fieldID:
      if field == DIM:
        dim = cs.unpack_int(DIM)
        if dim != 3: error("NWChem only performs 3d simulations")
      elif field == PERIODICITY:
        periodicity = cs.unpack(PERIODICITY,1)
        if basis_type == "ao":
          if periodicity[0] or periodicity[1] or periodicity[2]:
            error("NWChem AO basis wrapper only currently supports fully aperiodic systems")
        elif basis_type == "pw":
          if not periodicity[0] or not periodicity[1] or not periodicity[2]:
            error("NWChem PW basis wrapper only currently supports fully periodic systems")
      elif field == ORIGIN:
        origin = cs.unpack(ORIGIN,1)
      elif field == BOX:
        box_lmp = cs.unpack(BOX,1)
        if (basis_type == "pw"):
          if (box[0] != box_lmp[0] or box[4] != box_lmp[4] or box[8] != box_lmp[8]):
            error("NWChem wrapper mismatch in box dimensions")
      elif field == NATOMS:
        natoms_recv = cs.unpack_int(NATOMS)
        if natoms != natoms_recv:
          error("NWChem wrapper mismatch in number of atoms")
      elif field == COORDS:
        coords = cs.unpack(COORDS,1)

    if not origin or not box_lmp or not natoms or not coords:
      error("Required NWChem wrapper setup field not received");

  # STEP receive at each timestep of run or minimization
  # required fields: COORDS
  # optional fields: ORIGIN, BOX

  elif msgID == STEP:

    coords = []

    for field in fieldID:
      if field == COORDS:
        coords = cs.unpack(COORDS,1)

    if not coords: error("Required NWChem wrapper step field not received");

  else: error("NWChem wrapper received unrecognized message")

  # unpack coords from client
  # create NWChem input

  if basis_type == "ao":
    nwchem_input_write_ao(input_template,coords)
  elif basis_type == "pw":
    nwchem_input_write_pw(input_template,coords,box)

  # invoke NWChem

  subprocess.check_output(cp2kcmd,stderr=subprocess.STDOUT,shell=True)

  # process NWChem output
  
  i += 1
  log = "cp2k.out"

  if basis_type == "ao":
    energy,forces,charges = nwchem_read_ao(natoms,log)
    virial = [0,0,0,0,0,0]
  elif basis_type == "pw":
    energy,forces,virial,charges = nwchem_read_pw(log)

  #i += 1
  #archive = "./cp2k_logs/cp2k_" + str(i) + ".out"
  #shutil.copyfile(log,archive)
  #cmd_rm = "mv "+log+" ./cp2k_logs/cp2k_" + str(i) + ".out"
  #print(cmd_rm)
  #subprocess.check_output(cmd_rm,stderr=subprocess.STDOUT,shell=True)

  # return forces, energy, charges to client
  cs.send(msgID,4)
  cs.pack(FORCES,4,3*natoms,forces)
  cs.pack_double(ENERGY,energy)
  cs.pack(VIRIAL,4,6,virial)
  cs.pack(CHARGES,4,natoms,charges)

# final reply to client

cs.send(0,0)

# clean-up

del cs
