import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_T():
## Number of Directory and Directory Name ##
  ndir = 3;   dirnam = [[]]*ndir
  for n in range(ndir): dirnam[n] = './data%d' %(n+1)

## Computational Condition ##    
  dt = 50  # Time step in [a.u.]

####################################################################################################
# File to be read: md_eng.d
####################################################################################################

  dt *= 2.4188843265857e-5

## Array ##
  stp = np.array([]); temp = np.array([]);
   
## Read Data ##
  for n in range(ndir):
    filnam = '%s/md_eng.d' %(dirnam[n])
    fp = open(filnam, 'r')
    print('open: %s' %(filnam))
    while True:
      dat = fp.readline().split()   
      if len(dat) == 0: break
      if dat[0] == '#': continue
      stp = np.append(stp, float(dat[0]))
      temp = np.append(temp, float(dat[4]))
    fp.close()

## Console Output ##
  print('total number of data: %d' %(len(stp)))
  
## Plot Figure ##
  #plt.plot(stp*dt, temp)
  plt.plot(stp, temp)
  plt.title('Time Evolution of Temperature')
  #plt.xlabel('Time (ps)')
  plt.xlabel('MD step')
  plt.ylabel('Energy (eV)')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()
  
plot_T()
