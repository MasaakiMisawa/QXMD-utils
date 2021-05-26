import numpy as np
import numpy.linalg as nl

def get_trac():
## Number of Directory and Directory Name ##
  ndir = 11;   dirnam = [[]]*ndir
  dirnam[0] = '../data1'
  dirnam[1] = './data2'
  dirnam[2] = './data3'
  dirnam[3] = './data4'
  dirnam[4] = './data5'
  dirnam[5] = './data6'
  dirnam[6] = './data7'
  dirnam[7] = './data8'
  dirnam[8] = './data9'
  dirnam[9] = './dataA'
  dirnam[10] = './dataB'

  dirgrp = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # data group ( sum(dirgrp) = ndir )
  dirlbl = [i*0.05 for i in range(11)]        # list of first column ( len(dirlbl) = len(dirgrp) )

## Computational Conditions ##  
  avl = 100            # Average length (number of data)
  idn = 3              # Target plane ( 1: (100)  2: (010)  3: (001) )
  idm = 2              # Target direction ( 1: [100]  2: [010]  3: [001] )
  ofn = 'trac.dat'     # output file name

####################################################################################################
# Calculate average value of tractional vectors on the surfaces of simulation cell
#
# File to be read: (1) md_str.d  (2) md_cel.d
#
# Tractional vectors to be calculate:
#
#   t_n = (t_na, t_nb, t_nc)
#     where in a component of tractional vector t_nm (n, m = a, b, c),
#     n: define plane ( a -> (100), b -> (010), c -> (001) )
#     m: define direction (a -> vertical to (100), b -> vertical to (010), c -> vertical to (001))
#
#     if n != m, that means shear stress on the surface
#
# tractional vector is calculated based on the Cauchyi's law:
#
#  t = s n
#  t: traction vector
#  s: Cauchy stress tensor (symmetric 3x3 matrix)
#  n: normal vector parpendicular to a focused surface
#
####################################################################################################

  outn = dirgrp[0] - 1; cntg = 0; istr = [1,6,5,6,2,4,5,4,3]; fpo = open(ofn, 'w')
  t_a = np.array([[],[],[]]); t_b = np.array([[],[],[]]); t_c = np.array([[],[],[]])
  
  for n in range(ndir):
    filnam = '%s/md_cel.d' %(dirnam[n])
    fpc = open(filnam, 'r')
    fpc.readline(); fpc.readline()
    filnam = '%s/md_str.d' %(dirnam[n])
    fps = open(filnam, 'r')
    fps.readline(); fps.readline()
    
    while True:
      dac = np.array(fpc.readline().split(), dtype='float')   
      if len(dac) != 0:
        cel = np.array(dac[1:10]).reshape(3,3).T
        nx = np.matrix(np.cross(cel[:,1].T, cel[:,2].T) / nl.norm(np.cross(cel[:,1].T, cel[:,2].T)))
        ny = np.matrix(np.cross(cel[:,2].T, cel[:,0].T) / nl.norm(np.cross(cel[:,2].T, cel[:,0].T)))
        nz = np.matrix(np.cross(cel[:,0].T, cel[:,1].T) / nl.norm(np.cross(cel[:,0].T, cel[:,1].T)))
      das = np.array(fps.readline().split(), dtype='float')   
      if len(das) == 0: break
      str = np.array([])
      for i in istr: str = np.append(str, das[i])
      str = np.matrix(str.reshape(3,3))
      t_a = np.c_[t_a, str*nx.T]
      t_b = np.c_[t_b, str*ny.T]
      t_c = np.c_[t_c, str*nz.T]
    fpc.close()
    fps.close()
     
    if(n == outn):
      if idn == 1:
        print('%4.2f %9.6f' %(dirlbl[cntg], np.average(t_a[idm-1, t_a.shape[1] - avl : t_a.shape[1] ])))
        fpo.write('%4.2f %9.6f\n' %(dirlbl[cntg], np.average(t_a[idm-1, t_a.shape[1] - avl : t_a.shape[1] ])))
      if idn == 2:
        print('%4.2f %9.6f' %(dirlbl[cntg], np.average(t_b[idm-1, t_b.shape[1] - avl : t_b.shape[1] ])))
        fpo.write('%4.2f %9.6f\n' %(dirlbl[cntg], np.average(t_b[idm-1, t_b.shape[1] - avl : t_b.shape[1] ])))
      if idn == 3:
        print('%4.2f %9.6f' %(dirlbl[cntg], np.average(t_c[idm-1, t_c.shape[1] - avl : t_c.shape[1] ])))
        fpo.write('%4.2f %9.6f\n' %(dirlbl[cntg], np.average(t_c[idm-1, t_c.shape[1] - avl : t_c.shape[1] ])))
      cntg = cntg + 1
      if cntg < len(dirgrp): outn = outn + dirgrp[cntg]
      t_a = np.array([[],[],[]]); t_b = np.array([[],[],[]]); t_c = np.array([[],[],[]])
    
get_trac()
