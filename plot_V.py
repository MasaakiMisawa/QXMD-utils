import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_V():
## Number of Directory and Directory Name ##
  ndir = 3;   dirnam = [[]]*ndir
  for n in range(ndir): dirnam[n] = './data%d' %(n+1)

## Computational Conditions ##  
  nini = 500           # initial step number to calculate  
  dt = 50        # MD time step in [a.u.]

####################################################################################################

## Constant and Array ##
  ba = 0.529177210903; ba3 = ba*ba*ba; pi = np.arccos(-1)  
  dt *= 2.4188843265857e-5
  lini = 1; stp = np.array([]); vol = np.array([])
  a = np.array([]); b = np.array([]); c = np.array([])
  alp = np.array([]); bet = np.array([]); gam = np.array([])
  
## Read Data and Calculate Volume ##
  for n in range(ndir):
    filnam = '%s/md_cel.d' %(dirnam[n])
    fp = open(filnam, 'r')
    print('open: %s' %(filnam))
    while True:
      dat = fp.readline().split()   
      if dat == []: break
      if dat[0] == '#' or int(dat[0]) < nini: continue
      for i in range(10): dat[i] = float(dat[i])
      stp = np.append(stp, dat[0])
      cel = np.array(dat[1:10]).reshape(3,3).T
      a0 = np.linalg.norm(cel[:,0]); a = np.append(a, a0)
      b0 = np.linalg.norm(cel[:,1]); b = np.append(b, b0)
      c0 = np.linalg.norm(cel[:,2]); c = np.append(c, c0)
      alp = np.append(alp, np.arccos(np.dot(cel[:,1],cel[:,2])/(b0*c0))*180/pi)
      bet = np.append(bet, np.arccos(np.dot(cel[:,0],cel[:,2])/(a0*c0))*180/pi)
      gam = np.append(gam, np.arccos(np.dot(cel[:,1],cel[:,0])/(b0*a0))*180/pi)
      vol0 = np.dot(cel[:,0],np.cross(cel[:,1],cel[:,2])); vol = np.append(vol, vol0)
      if lini == 1: bvol = vol0; lini = 0 
    fp.close()
   
## Console Output ## 
  print('total number of data       : %d' %(len(stp)))
  print('average a (angstrom)       : %f' %(np.average(a)*ba))
  print('average b (angstrom)       : %f' %(np.average(b)*ba))
  print('average c (angstrom)       : %f' %(np.average(c)*ba))
  print('average alpha (degree)     : %f' %(np.average(alp)))
  print('average beta (degree)      : %f' %(np.average(bet)))
  print('average gamma (degree)     : %f' %(np.average(gam)))
  print('average volume (angstrom^3): %f' %(np.average(vol*ba3)))
  print('variance of volume         : %f' %(np.var(vol*ba3)))
  print('variance/average of volume : %f' %(np.var(vol*ba3)/np.average(vol*ba3)))
  
## Output File ##
  fp = open('Volume.dat', 'w')
  fp.write('# step, volume (angstrom^3)\n')
  for i in range(len(vol)): fp.write('%5d %13.6f\n' %(stp[i],vol[i]*ba3))
  
## Plot Figure ##
  plt.title('Time Evolution of Volume')
  #plt.plot(stp*dt, vol*ba3*1.e-3)
  #plt.xlabel('Time (ps)'); plt.ylabel('Volume (nm$^{3}$)')
  plt.plot(stp, vol*ba3*1.e-3)  
  plt.xlabel('MD step'); plt.ylabel('Volume (nm$^{3}$)')
  #plt.plot(stp, vol/bvol)  
  #plt.xlabel('MD step'); plt.ylabel('Specific volume')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()

plot_V()
