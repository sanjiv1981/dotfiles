#!/usr/bin/python
from sys import *

if len(argv) == 1:
  print "Syntax: numelectrons.py <FILES>"
  print "This script is designed to calculate the number of electrons in the first"
  print " frame of an XYZ file"
  exit(1)


for ifile in argv[1:]:
  
  Nelectron=0
  with open(ifile,'r') as f:
    Natom=int(f.readline())
    f.readline()
    for i in range(Natom):
      ChemSym = f.readline().split()[0]
      
      if   ChemSym == 'H':
        Nelectron += 1
      elif ChemSym == 'He':
        Nelectron += 2
      elif ChemSym == 'Li':
        Nelectron += 3
      elif ChemSym == 'Be':
        Nelectron += 4
      elif ChemSym == 'B':
        Nelectron += 5
      elif ChemSym == 'C':
        Nelectron += 6
      elif ChemSym == 'N':
        Nelectron += 7
      elif ChemSym == 'O':
        Nelectron += 8
      elif ChemSym == 'F':
        Nelectron += 9
      elif ChemSym == 'Ne':
        Nelectron += 10
      elif ChemSym == 'Na':
        Nelectron += 11
      elif ChemSym == 'Mg':
        Nelectron += 12
      elif ChemSym == 'Al':
        Nelectron += 13
      elif ChemSym == 'Si':
        Nelectron += 14
      elif ChemSym == 'P':
        Nelectron += 15
      elif ChemSym == 'S':
        Nelectron += 16
      elif ChemSym == 'Cl':
        Nelectron += 17
      elif ChemSym == 'Ar':
        Nelectron += 18
      else:
        print 'Chem Symbol ',ChemSym, ' not defined'
        raise

  print Nelectron,ifile
