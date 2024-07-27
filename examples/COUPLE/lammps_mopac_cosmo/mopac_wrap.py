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
#vaspcmd = "mopac in.mop"
vaspcmd = "/mnt/d/mopac-22.1.1/bin/mopac in.mop"

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

def vasp_setup(poscar): # poscar = mop format

  ps = open(poscar,'r').readlines()
  
  box = []
  types_symbol = []
  QMatoms = 0
  ntypes = 0
  for line in ps[2:]:
    words = line.split()
    if len(words) >= 1:
      if words[0] == "Tv":
        box += [words[1],words[2],words[3]]
      else:
        types_symbol += [words[0]]
        QMatoms += 1
  
  #print(types_symbol)
  ntypes = len(types_symbol)
  
  if len(box) == 0:
    box = [0.0,0.0,0.0,
           0.0,0.0,0.0,
           0.0,0.0,0.0]

  return QMatoms,ntypes,box,types_symbol
  
# -------------------------------------
# write a new POSCAR file for VASP

def poscar_write(poscar,QMatoms,ntypes,types,coords,box,types_symbol):

  psold = open(poscar,'r').readlines()
  psnew = open("in.mop",'w')
  
  words = psold[0].split()
  sttype = words[1]

  # header (keywords and comment line)
  
  psnew.write(psold[0]) # comment line
  psnew.write(psold[1]) # comment line
  psnew.write(psold[2]) # comment line
  
  # per-atom coords
  # grouped by types
  
  t = 0
  for i in range(QMatoms):
    t = types_symbol[i]
    x = coords[3*i+0]
    y = coords[3*i+1]
    z = coords[3*i+2]
    aline = "%s  %g 1 %g 1 %g 1 \n" % (t,x,y,z)
    psnew.write(aline)

  #if box[0] > 0.0 or box[1] > 0.0 or box[2] > 0.0: psnew.write("Tv %g %g %g \n" % (box[0], box[1], box[2])) # a
  #if box[3] > 0.0 or box[4] > 0.0 or box[5] > 0.0: psnew.write("Tv %g %g %g \n" % (box[3], box[4], box[5])) # b
  #if box[6] > 0.0 or box[7] > 0.0 or box[8] > 0.0: psnew.write("Tv %g %g %g \n" % (box[6], box[7], box[8])) # c

  psnew.close()
  
# -------------------------------------
# write a new POSCAR file for VASP

#def charge_write(QMatoms,natoms,ntypes,types,coords,lmp_charges,types_symbol):
#
#  psnew = open("mol.in",'w')
#  
#  psnew.write("empty_line \n")
#  aline = "%d %d     # of qmmm atoms, # of link atoms in Region I \n" % (QMatoms,0)
#  psnew.write(aline)
#  
#  #coulomb_constants = (1.602176634e-19)**2/(4*3.1415926*8.85418781e-12)/(1e-10) # Nm/q^2 = J/q^2
#  #coulomb_constants = 2.30662e-18 / 1.60217663E-19 # J/q^2 to eV/q^2
#  #coulomb_constants = 14.39677 * 23.062063375661 # eV/q^2 to kcal/mol/q^2
#  coulomb_constants = 332.0193 # kcal/mol/q^2
#  phi = []
#  for i in range(QMatoms): # QM
#    t = types_symbol[i]
#    x = coords[3*i+0]
#    y = coords[3*i+1]
#    z = coords[3*i+2]
#    q = lmp_charges[i]
#    #---------------------------------------------
#    coule = 0
#    for j in range(QMatoms, natoms): # MM
#      xj = coords[3*j+0]
#      yj = coords[3*j+1]
#      zj = coords[3*j+2]
#      qj = lmp_charges[j]
#      r = ((xj-x)**2+(yj-y)**2+(zj-z)**2)**0.5
#      coule += qj/r
#    #---------------------------------------------
#    aline = "%s  %g %g %g  %g\n" % (t,x,y,z,coulomb_constants*coule)
#    psnew.write(aline)
#
#  psnew.close()

# -------------------------------------
# read a VASP output vasprun.xml file
# uses ElementTree module
# see https://docs.python.org/2/library/xml.etree.elementtree.html

def vasprun_read():
  
  #-------------------------------------------------------------------
  # detailed.out version
  #-------------------------------------------------------------------

  qout = []
  fout = []
  sout = []
  stensor = []

  eunitconv=0.159466838598749E-02/0.367493245336341e-01 # kcal/mol to eV
  funitconv=0.159466838598749E-02/0.367493245336341e-01 # kcal/mol/Angstrom to eV/Angstrom
  
  with open('in.aux', 'r') as file:
    #--------------------------------------------------------------
    # HEAT_OF_FORMATION
    start_reading = False
    for line in file:
      line_list = re.split(r"[:=\[\]]", line.strip("\n").strip())
      if 'HEAT_OF_FORMATION' in line_list[0].strip():
        eout = float(line_list[2].replace("D", "E"))*eunitconv
        break
    #print(eout)
    #--------------------------------------------------------------
    # ATOM_CHARGES
    start_reading = False
    for line in file:
      line_list = re.split(r"[:=\[\]]", line.strip("\n").strip())
      if start_reading and line_list[0] == 'AO_CHARGES':
        break
      if start_reading:
        qxyz = line.strip().split()
        qxyz = [float(value) for value in qxyz]
        qout += qxyz
      if 'ATOM_CHARGES' in line_list[0]:
        start_reading = True
        continue
    #print(qout)
    #--------------------------------------------------------------
    # GRADIENTS
    start_reading = False
    for line in file:
      line_list = re.split(r"[:=\[\]]", line.strip("\n").strip())
      if start_reading and line_list[0] == 'OVERLAP_MATRIX':
        break
      if start_reading:
        fxyz = line.strip().split()
        fxyz = [float(value)*funitconv for value in fxyz]
        fout += fxyz
      #else:
      if 'GRADIENTS' in line_list[0]:
        start_reading = True
        continue
    #print(fout)
    #--------------------------------------------------------------

#  sunitconv=1.0/0.339893208050290E-13
#  start_reading = False
#  with open('results.tag', 'r') as file:
#    for line in file:
#      line_list = re.split(r"[:=\[\]]", line.strip("\n").strip())
#      if start_reading and line_list[0] == 'forces_ext_charges':
#        break
#      if start_reading and line_list[0] == 'cell_volume':
#        break
#      if start_reading and line_list[0] == 'orbital_charges':
#        break
#      if start_reading:
#        sxyz = line.strip().split()
#        sxyz = [float(value)*sunitconv for value in sxyz]
#        stensor.append(sxyz)
#      if 'stress' in line_list[0]:
#        start_reading = True
#        continue
#  if start_reading:
#    sxx = stensor[0][0]
#    syy = stensor[1][1]
#    szz = stensor[2][2]
#    # symmetrize off-diagonal components
#    sxy = 0.5 * (stensor[0][1] + stensor[1][0])
#    sxz = 0.5 * (stensor[0][2] + stensor[2][0])
#    syz = 0.5 * (stensor[1][2] + stensor[2][1])
#    sout = [sxx,syy,szz,sxy,sxz,syz]
#  else:
#    sout = [0.0,0.0,0.0,0.0,0.0,0.0]
  #print(sout)
  sout = [0.0,0.0,0.0,0.0,0.0,0.0]
  #-------------------------------------------------------------------
  
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
QMatoms,ntypes,box,types_symbol = vasp_setup(poscar_template)

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
  
  poscar_write(poscar_template,QMatoms,ntypes,types,coords,box,types_symbol)
  
  # create charges.dat file
  #charge_write(QMatoms,natoms,ntypes,types,coords,lmp_charges,types_symbol)

  # invoke DFTB+
  
  #print("\nLaunching DFTB+ ...")
  #print(vaspcmd)
  subprocess.check_output(vaspcmd,stderr=subprocess.STDOUT,shell=True)
  
  # process VASP output
    
  energy,forces,virial,charges = vasprun_read()
  
  # coulomb interaction (forces and energies) of MM area from QM area
  #coulomb_constants = (1.602176634e-19)**2/(4*3.1415926*8.85418781e-12)/(1e-10) # Nm/q^2 = J/q^2
  #coulomb_constants = 2.30662e-18 / 1.60217663E-19 # J/q^2 to eV/q^2
  #coulomb_constants = 14.39677 * 23.062063375661 # eV/q^2 to kcal/mol/q^2
  #coulomb_constants = 332.0193 # kcal/mol/q^2
  #funitconv=0.159466838598749E-02/0.367493245336341e-01 # kcal/mol/Angstrom to eV/Angstrom
  #for j in range(QMatoms, natoms): # MM
  #  xj = coords[3*j+0]
  #  yj = coords[3*j+1]
  #  zj = coords[3*j+2]
  #  qj = lmp_charges[j]
  #  #-----------------------------------------------------
  #  fcq = 0.0
  #  fxj = 0.0
  #  fyj = 0.0
  #  fzj = 0.0
  #  for i in range(QMatoms): # QM
  #    xi = coords[3*i+0]
  #    yi = coords[3*i+1]
  #    zi = coords[3*i+2]
  #    qi = lmp_charges[i]
  #    r = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**0.5
  #    fcq = funitconv*(coulomb_constants*qi*qj/(r**2))
  #    fxj += -fcq*( (xi-xj)/r ) # eV/Angstrom
  #    fyj += -fcq*( (yi-yj)/r ) # eV/Angstrom
  #    fzj += -fcq*( (zi-zj)/r ) # eV/Angstrom
  #    #energy += funitconv*(coulomb_constants*qi*qj/r) # eV
  #    # the electron-electron part of Hartree-Fock equation
  #    # for i in range(nq):
  #    #   for j in range(i+1,nq):
  #    #     energy += funitconv*(coulomb_constants*qi*qj/r)
  #    # This part is already included in the form of potential in 
  #    # the last column of the "mol.in" file in MOPAC, 
  #    # so no calculation is necessary.
  #    # Please correct me if I'm wrong.
  #  #-----------------------------------------------------
  #  forces += [fxj,fyj,fzj]
  
  #for i in range(QMatoms,natoms):
  #  charges += [float(lmp_charges[i])]
  
  # debag
  #print(forces)
  #print(virial)
  #print(charges)

  # convert VASP kilobars to bars

  #for i,value in enumerate(virial): virial[i] *= 1.0/(100.0*1000.0) # Pa to bars (1 bar = 0.1 MPa = 100 kPa)
  
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
