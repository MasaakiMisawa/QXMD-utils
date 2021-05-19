import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_HUG():
## Number of Directory and Directory Name ##
  ndir = 1;   dirnam = [[]]*ndir
  for n in range(ndir): dirnam[n] = './MSST/data2'

## Computational Conditions ##  
  nini = 500           # Initial step 
  nend = 1e9           # Final step
  nave = 1000          # Average range [step]
  v0 = 0               # Initial volume (< 0: volume in bohr^3,  >=0: referance step)
  t0 = -1              # Time of the first plot (-1: depend on step number)
  dt = 0.000242        # MD time step in [a.u.]

####################################################################################################

## Constant and Array ##
  pi = np.arccos(-1);  dt *= 2.4188843265857e-5
  lini = 1; bt = 0; bvol = 1; stp = np.array([]) 
  vol = np.array([]); vel = np.array([]); str = np.array([])
  rvol = np.array([]); rvel = np.array([]); rstr = np.array([])
  
## Read Data and Calculate Volume ##
  for n in range(ndir):
    fpc = open('%s/md_cel.d' %(dirnam[n]), 'r')
    for i in range(2): fpc.readline().split() 
    fph = open('%s/md_hug.d' %(dirnam[n]), 'r')
    for i in range(4): fph.readline().split() 
    while True:
      dtc = np.array(fpc.readline().split(), dtype='float')   
      if len(dtc) == 0: break
      dth = np.array(fph.readline().split(), dtype='float')
      if int(dth[0]) < nini: continue
      if lini == 1 and t0 >= 0:
        bt = dtc[0]; lini = 0
      cel = np.array(dtc[1:10]).reshape(3,3).T
      vol0 = np.dot(cel[:,0],np.cross(cel[:,1],cel[:,2])); vol = np.append(vol, vol0)
      if v0 == int(dtc[0]): 
        bvol = vol0
      stp = np.append(stp, dth[0]); 
      vel = np.append(vel, dth[1])
      str = np.append(str, dth[2])
    fpc.close()
    fph.close()
  
  if len(stp) - nave + 1 <= 0:
    ptint('Error: too small nave'); return 1
  vel *= 1.e-3  # m/s -> km/s
  for i in range(len(stp) - nave + 1):
    avol = np.average(vol[i:i+nave]); rvol = np.append(rvol, avol)
    avel = np.average(vel[i:i+nave]); rvel = np.append(rvel, avel)
    astr = np.average(str[i:i+nave]); rstr = np.append(rstr, astr)
    
## Console Output ## 
  print('Average value for latest %d step' %(nave))
  print('V/V0: %f' %(avol/bvol))
  print('Up  : %f' %(avel))
  print('P   : %f' %(astr))
  
## Plot Figure ##
  plt.title('Time Evolution of Hugoniot Pressure')
  plt.plot((stp - bt)*dt, str)  
  plt.plot((stp[len(stp)-nave-1:] - bt)*dt, rstr)  
  plt.xlabel('Time (ps)'); plt.ylabel('Pressure (GPa)')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()

  plt.title('Time Evolution of Particle Velocity')
  plt.plot((stp - bt)*dt, vel)  
  plt.plot((stp[len(stp)-nave-1:] - bt)*dt, rvel)  
  plt.xlabel('Time (ps)'); plt.ylabel('Particle velocity (km/s)')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()

  plt.title('Time Evolution of Specific Volume')
  plt.plot((stp - bt)*dt, vol/bvol)  
  plt.plot((stp[len(stp)-nave-1:] - bt)*dt, rvol/bvol)  
  plt.xlabel('Time (ps)'); plt.ylabel('Specific Volume')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()

plot_HUG()
