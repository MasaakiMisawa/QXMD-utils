import numpy as np
#mpl.use('SVG')
import matplotlib.pyplot as plt

def plot_E():
## Number of Directory and Directory Name ##
  ndir = 3;   dirnam = [[]]*ndir
  for n in range(ndir): dirnam[n] = './data%d' %(n+1)
  
## Computational Condition ##    
  dt = 50    # Time step in [a.u.]
  
####################################################################################################
# File to be read: md_eng.d
####################################################################################################

## Constant and Array ##  
  hev = 27.211396127707; stp = np.array([]); 
  dt *= 2.4188843265857e-5
  etot = np.array([]); epot = np.array([]); 
  base = np.array([0.,0.,0.]); lini = 1
  
## Read Data ##
  for n in range(ndir):
    filnam = '%s/md_eng.d' %(dirnam[n])
    fp = open(filnam, 'r')
    print('open: %s' %(filnam))
    while True:
      dat = fp.readline().split()   
      if dat == []: break
      if dat[0] == '#': continue
      if lini == 1:
        for i in range(1, 3): base[i] = float(dat[i])
        lini = 0
      stp = np.append(stp, float(dat[0]))
      etot = np.append(etot, float(dat[1]))
      epot = np.append(epot, float(dat[2]))
    fp.close()
    
## Console Output ##  
  print('total number of data: %d' %(len(stp)))
  
## Plot Figure ##
  #plt.plot(stp*dt, (etot-base[1])*hev, label = 'Total energy')
  #plt.plot(stp*dt, (epot-base[2])*hev, label = 'Potential energy')
  #plt.plot(stp, (etot-base[1])*hev, label = 'Total energy')
  #plt.plot(stp, (epot-base[2])*hev, label = 'Potential energy')
  plt.plot(stp, etot*hev, label = 'Total energy')
  plt.plot(stp, epot*hev, label = 'Potential energy')
  plt.legend()
  plt.title('Time Evolution of Energy')
  #plt.xlabel('Time (ps)')
  plt.xlabel('MD step')
  plt.ylabel('Energy (eV)')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()
  #plt.savefig('test.svg')
  
plot_E()
