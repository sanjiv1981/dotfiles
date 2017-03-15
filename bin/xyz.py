#################################################################################
# A Set of utilities for reading, writing, and manipulating atomic configurations
#################################################################################
#
#  This set of utilities uses the XYZ file format.
#  The top level routines available are:
#
# DensifyChain(inputchain,Nnew=1):
# readXYZ(filename):
# writeXYZ(config,filename):
# writechain(chain, filename):
# readchain(filename,Nimage):
# createchain(start,final,Nimage):
# MoveToOrigin(config,atomID):
# RotateAxisAngle(config, unit_axis, angle):
# RotateVectorToDirection(config, vector, direction):
# reflectXZ(config):
# reflectXY(config):
# reflectYZ(config):
# rotateX(config, angle):
# rotateY(config, angle):
# rotateZ(config, angle):


import numpy
import sys
import copy

#==================================================================================
class ConfigClass:
#   Natom    : Scalar,      Number of atoms
#   E        : Scalar,      Energy of configuration. Initialized to zero
#   Gradient : Numpy array, Gradient, initialized to zero
#   ChemSym  : List,        Chemical symbol of each atom
#   Coord    : Numpy array, XYZ coordinates of each atom. First index is atom number,
#                           Second index is coordinate direction
#==================================================================================
  def __init__(self):
    self.ChemSym  = []
    self.Coord    = numpy.array([]).reshape(0,3)
    self.Natom    = 0
    self.E        = 0.0
    self.Gradient = 0*self.Coord

#========================================================================================================  
def DensifyChain(inputchain,Nnew=1):
# Routine that densifies a chain, by inserting Nnew images
# in between every successive pair of images of the original
# chain. The end points of the chain are kept fixed
#========================================================================================================  
  Nimage = len(inputchain)
  print Nimage
  
  Nnewimage = (Nimage-1)*(Nnew+1)+1
  newchain = [None]*Nnewimage
  
  k = 0
  for i in range(Nimage-1):
    dr = (inputchain[i+1].Coord - inputchain[i].Coord)/(Nnew+1)
    for j in range(Nnew+1):
      newchain[k] = copy.deepcopy(inputchain[i])
      newchain[k].Coord = inputchain[i].Coord + j*dr
      k+=1
  newchain[k] = copy.deepcopy(inputchain[-1])
  return(newchain)
  #sys.exit()
  
#========================================================================================================  
def readXYZ(filename):
#  Routine that reads in XYZ file format.
#  Expected format is
#     First line         : Number of atoms
#     Second line        : Comment
#     Third line onwards : One line per atom, ChemSym X Y Z
#========================================================================================================  

  try:
    ifile=open(filename,'r')
  except:
    print '------> Cannot open ',filename,' <------'
    raise
  config = ConfigClass()
  
  config.Natom = int(ifile.readline().strip())
  sdum  = ifile.readline()
  for iatom in range(config.Natom):
    line = ifile.readline().split()
    config.ChemSym.append(line[0])
    config.Coord = numpy.append(config.Coord, [numpy.array(map(float,line[1:]))], axis=0)
  ifile.close()
  config.Gradient = 0*config.Coord
  config.E = 0
  return(config)

#========================================================================================================  
def writeXYZ(config,filename):
#  Routine that writes a config in XYZ file format.
#  Output format is
#     First line         : Number of atoms
#     Second line        : Comment
#     Third line onwards : One line per atom, ChemSym X Y Z
#========================================================================================================  
  try:
    ofile = open(filename,'w')
  except:
    print '------> Cannot open file ',filename,' for writing <------'
    raise
  
  ofile.write('%d\n' %(config.Natom))
  ofile.write('XYZ format file, created by XYZ Utility\n')
  for iatom in range(config.Natom):
    ofile.write(config.ChemSym[iatom])
    for j in range(3):
      ofile.write('%15.8f ' %(config.Coord[iatom][j]))
    ofile.write('\n')
  ofile.close()
  
  
#========================================================================================================  
def writechain(chain, filename):
#========================================================================================================  
  try:
    ofile = open(filename,'w')
  except:
    print '------> Cannot open chain file ',filename,' for writing <------'
    raise
  
  Nimage = len(chain)
  for iimg in range(Nimage):
    ofile.write('%d\n' %(chain[iimg].Natom))
    ofile.write('Image Number %d\n' %(iimg))
    for iatom in range(chain[iimg].Natom):
      ofile.write(chain[iimg].ChemSym[iatom]),
      for j in range(3):
	ofile.write("%15.8f " %(chain[iimg].Coord[iatom][j]))
      ofile.write("\n")
  ofile.close()
  
#========================================================================================================  
def readchain(filename,Nimage):
#========================================================================================================  
  try:
    ifile  = open(filename,'r')
  except:
    print '------> Cannot open chain file ',filename,' for reading <------'
    raise
  
  returnconfig = [None]*Nimage
  for iimg in range(Nimage):
    returnconfig[iimg] = ConfigClass()
    Natom = int(ifile.readline())
    
    ldum = ifile.readline()
    for iatom in range(Natom):
      line = ifile.readline().split()
      returnconfig[iimg].ChemSym.append(line[0])
      returnconfig[iimg].Coord = numpy.append(returnconfig[iimg].Coord,[numpy.array(map(float,line[1:]))],axis=0)
    returnconfig[iimg].Gradient = 0*returnconfig[iimg].Coord
    returnconfig[iimg].Natom    = Natom
    returnconfig[iimg].E        = 0
  ifile.close()
  return(returnconfig)

#========================================================================================================  
def createchain(start,final,Nimage):
#
#  Create a chain-of-states, by linearly interpolating between start and final configs
#========================================================================================================  
  Natom = start.Natom
  if Natom != final.Natom:
    print '------> start and final configs do not have the same number of atoms <------'
    raise
  
  chain = [None]*Nimage
  dr = (final.Coord - start.Coord)/(Nimage-1)
  for iimg in range(Nimage):
    chain[iimg] = ConfigClass()
    chain[iimg].ChemSym = start.ChemSym
    chain[iimg].Coord   = start.Coord + iimg*dr
    chain[iimg].Natom   = Natom
  return(chain)

#========================================================================================================  
def MoveToOrigin(config,atomID):
#
#  Translate the entire configuration, so that <atomID> is at the origin.
#========================================================================================================  
  Natom    = config.Natom
  shiftvec = copy.deepcopy(config.Coord[atomID])
  
  for iatom in range(Natom):
    config.Coord[iatom] -= shiftvec
      
#========================================================================================================  
def RotateAxisAngle(config, unit_axis, angle):
#
#  Rotate the configuration about the origin, about <unit_axis>, by <angle> (radians)
#========================================================================================================  
  Natom = config.Natom
  u = unit_axis[0]
  v = unit_axis[1]
  w = unit_axis[2]
  ct = numpy.cos(angle)
  st = numpy.sin(angle)
  
  for iatom in range(Natom):
    x = copy.deepcopy(config.Coord[iatom][0])
    y = copy.deepcopy(config.Coord[iatom][1])
    z = copy.deepcopy(config.Coord[iatom][2])
    config.Coord[iatom][0] = u*(u*x+v*y+w*z)*(1-ct) + x*ct + (-w*y+v*z)*st
    config.Coord[iatom][1] = v*(u*x+v*y+w*z)*(1-ct) + y*ct + ( w*x-u*z)*st
    config.Coord[iatom][2] = w*(u*x+v*y+w*z)*(1-ct) + z*ct + (-v*x+u*y)*st
    
#========================================================================================================  
def RotateVectorToDirection(config, vector, direction):
#
#  Rotate the configuration, so that <vector> is pointed along <direction>
#========================================================================================================  
  
  # Normalize input vectors
  unit_vector    = vector    / numpy.sqrt(numpy.dot(vector,vector))
  unit_direction = direction / numpy.sqrt(numpy.dot(direction,direction))
  
  # Find unit axis of rotation
  axis = numpy.cross(unit_vector, unit_direction)
  unit_axis = axis / numpy.sqrt(numpy.dot(axis,axis))
  
  # Find angle
  theta = numpy.arccos(numpy.dot(unit_vector,unit_direction) )
  
  RotateAxisAngle(config,unit_axis,theta)
  
#========================================================================================================  
def reflectXZ(config):
#
#  Preform a reflection through the XZ plane
#========================================================================================================  
  Natom = config.Natom
  for i in range(Natom):
    config.Coord[i][1] = -config.Coord[i][1]

#========================================================================================================  
def reflectXY(config):
#  Preform a reflection through the XY plane
#========================================================================================================  
  Natom = config.Natom
  for i in range(Natom):
    config.Coord[i][2] = -config.Coord[i][2]
    
#========================================================================================================  
def reflectYZ(config):
#  Preform a reflection through the YZ plane
#========================================================================================================  
  Natom = config.Natom
  for i in range(Natom):
    config.Coord[i][0] = -config.Coord[i][0]

#========================================================================================================  
def rotateX(config, angle):
#========================================================================================================  
  RotateAxisAngle(config,numpy.array([1,0,0],dtype=float),angle)

#========================================================================================================  
def rotateY(config, angle):
#========================================================================================================  
  RotateAxisAngle(config,numpy.array([0,1,0],dtype=float),angle)

#========================================================================================================  
def rotateZ(config, angle):
#========================================================================================================  
  RotateAxisAngle(config,numpy.array([0,0,1],dtype=float),angle)
