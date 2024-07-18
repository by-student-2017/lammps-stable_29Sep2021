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
vaspcmd = "dftb+"

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

def vasp_setup(poscar):

  ps = open(poscar,'r').readlines()

  # box size

  words = ps[2].split()
  xbox = float(words[0])
  words = ps[3].split()
  ybox = float(words[1])
  words = ps[4].split()
  zbox = float(words[2])
  box = [xbox,ybox,zbox]

  ntypes = 0
  natoms = 0
  words = ps[6].split()
  for word in words:
    if word == '#': break
    ntypes += 1
    natoms += int(word)
  
  return natoms,ntypes,box
  
# -------------------------------------
# write a new POSCAR file for VASP

def poscar_write(poscar,natoms,ntypes,types,coords,box):

  psold = open(poscar,'r').readlines()
  psnew = open("POSCAR",'w')

  # header, including box size
  
  psnew.write(psold[0])
  psnew.write(psold[1])
  psnew.write("%g %g %g\n" % (box[0],box[1],box[2]))
  psnew.write("%g %g %g\n" % (box[3],box[4],box[5]))
  psnew.write("%g %g %g\n" % (box[6],box[7],box[8]))
  psnew.write(psold[5])
  psnew.write(psold[6])

  # per-atom coords
  # grouped by types

  psnew.write("Cartesian\n")

  for itype in range(1,ntypes+1):
    for i in range(natoms):
      if types[i] != itype: continue
      x = coords[3*i+0]
      y = coords[3*i+1]
      z = coords[3*i+2]
      aline = "  %g %g %g\n" % (x,y,z)
      psnew.write(aline)

  psnew.close()

# -------------------------------------
# read a VASP output vasprun.xml file
# uses ElementTree module
# see https://docs.python.org/2/library/xml.etree.elementtree.html

def vasprun_read():
  
  #-------------------------------------------------------------------
  # detailed.out version
  #-------------------------------------------------------------------

  start_reading = False
  with open('results.tag', 'r') as file:
    for line in file:
      flag_line = line.split(":")
      flag_line[0] = flag_line[0].strip()
      if 'mermin_energy' in flag_line[0]:
        start_reading = True
        continue
      end_text = line.strip().split(":")
      if start_reading and flag_line[0] == 'extrapolated0_energy':
        break
      if start_reading:
        values = line.strip().split()
        eout = float(values[0])
  #print(eout)

  fout = []
  sout = []
  stensor = []
  qout = []
  
  start_reading = False
  with open('results.tag', 'r') as file:
    for line in file:
      flag_line = line.split(":")
      flag_line[0] = flag_line[0].strip()
      if 'forces' in flag_line[0]:
        start_reading = True
        continue
      if start_reading and flag_line[0] == 'stress':
        break
      if start_reading:
        fxyz = line.strip().split()
        fxyz = [float(value) for value in fxyz]
        fout += fxyz
  #print(fout)

  start_reading = False
  with open('results.tag', 'r') as file:
    for line in file:
      flag_line = line.split(":")
      flag_line[0] = flag_line[0].strip()
      if 'stress' in flag_line[0]:
        start_reading = True
        continue
      if start_reading and flag_line[0] == 'cell_volume':
        break
      if start_reading:
        sxyz = line.strip().split()
        sxyz = [float(value) for value in sxyz]
        stensor.append(sxyz)
  sxx = stensor[0][0]
  syy = stensor[1][1]
  szz = stensor[2][2]
  # symmetrize off-diagonal components
  sxy = 0.5 * (stensor[0][1] + stensor[1][0])
  sxz = 0.5 * (stensor[0][2] + stensor[2][0])
  syz = 0.5 * (stensor[1][2] + stensor[2][1])
  sout = [sxx,syy,szz,sxy,sxz,syz]
  #print(sout)

  start_reading = False
  with open('results.tag', 'r') as file:
    for line in file:
      flag_line = line.split(":")
      flag_line[0] = flag_line[0].strip()
      if 'gross_atomic_charges' in flag_line[0]:
        start_reading = True
        continue
      if start_reading and flag_line[0] == 'dipole_moments':
        break
      if start_reading:
        qxyz = line.strip().split()
        qxyz = [float(value) for value in qxyz]
        qout += qxyz
  #print(qout)
  #-------------------------------------------------------------------
  
  return eout,fout,sout,qout
  
  #-------------------------------------------------------------------
  # detailed.out version
  #-------------------------------------------------------------------
#  search_term="Total Mermin free energy:"
#  with open('detailed.out', 'r') as file:
#    for line in file:
#      if re.search(search_term, line):
#         values = line.strip().split()
#         eout = float(values[6])
  #print(eout)
#
#  fout = []
#  sout = []
#  stensor = []
#  
#  start_reading = False
#  with open('detailed.out', 'r') as file:
#    for line in file:
#      if 'Total Forces' in line:
#        start_reading = True
#        continue
#      if start_reading and line.strip() == '':
#        break
#      if start_reading:
#        force = line.strip().split()
#        fxyz = force[1:]
#        fxyz = [float(value) for value in fxyz]
#        fout += fxyz
  #print(fout)
#
#  start_reading = False
#  with open('detailed.out', 'r') as file:
#    for line in file:
#      if 'Total stress tensor' in line:
#        start_reading = True
#        continue
#      if start_reading and line.strip() == '':
#        break
#      if start_reading:
#        sxyz = line.strip().split()
#        sxyz = [float(value) for value in sxyz]
#        stensor.append(sxyz)
#  sxx = stensor[0][0]
#  syy = stensor[1][1]
#  szz = stensor[2][2]
  # symmetrize off-diagonal components
#  sxy = 0.5 * (stensor[0][1] + stensor[1][0])
#  sxz = 0.5 * (stensor[0][2] + stensor[2][0])
#  syz = 0.5 * (stensor[1][2] + stensor[2][1])
#  sout = [sxx,syy,szz,sxy,sxz,syz]
  #print(sout)
  #-------------------------------------------------------------------
  
#  return eout,fout,sout

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

natoms,ntypes,box = vasp_setup(poscar_template)

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
    
    for field in fieldID:
      if field == DIM:
        dim = cs.unpack_int(DIM)
        if dim != 3: error("VASP only performs 3d simulations")
      elif field == PERIODICITY:
        periodicity = cs.unpack(PERIODICITY,1)
        if not periodicity[0] or not periodicity[1] or not periodicity[2]:
          error("VASP wrapper only currently supports fully periodic systems")
      elif field == ORIGIN:
        origin = cs.unpack(ORIGIN,1)
      elif field == BOX:
        box = cs.unpack(BOX,1)
      elif field == NATOMS:
        natoms_recv = cs.unpack_int(NATOMS)
        if natoms != natoms_recv:
          error("VASP wrapper mis-match in number of atoms")
      elif field == NTYPES:
        ntypes_recv = cs.unpack_int(NTYPES)
        if ntypes != ntypes_recv:
          error("VASP wrapper mis-match in number of atom types")
      elif field == TYPES:
        types = cs.unpack(TYPES,1)
      elif field == COORDS:
        coords = cs.unpack(COORDS,1)

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
  
  poscar_write(poscar_template,natoms,ntypes,types,coords,box)

  # invoke DFTB+
  
  #print("\nLaunching DFTB+ ...")
  #print(vaspcmd)
  subprocess.check_output(vaspcmd,stderr=subprocess.STDOUT,shell=True)
  
  # process VASP output
    
  energy,forces,virial,charges = vasprun_read()

  # convert VASP kilobars to bars

  for i,value in enumerate(virial): virial[i] *= 1000.0
    
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
