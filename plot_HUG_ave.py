import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_avHUG():
## Number of Directory and Directory Name ##
  ndir = 16;   dirnam = [[]]*ndir
  for n in range(ndir): dirnam[n] = './MSST/%d' %(6000+100*n) 

## Computational Conditions ##  
  nini = 500           # Initial step 
  nend = 1e9           # Final step
  nave = 1500          # Average range [step]
  v0   = -23273.904371 # Initial volume (< 0: volume in bohr^3,  >=0: referance step)
  t0   = -1            # Time of the first plot (-1: depend on step number)
  ldmp = True         # Make datafile

####################################################################################################
# Files to be read: (1) md_cel.d (2) md_hug.d
####################################################################################################

## Constant and Array ##
  lini = 1; bt = 0; bvol = -v0;  
  vol_sd = np.array([]); vel_sd = np.array([]); str_sd = np.array([])
  avol = np.array([]); avel = np.array([]); astr = np.array([])
  vshk = np.array([6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5])

## Reference Data ## 
  vshk_ref = np.array([6.000, 6.100, 6.200, 6.500, 6.600, 6.700, 7.000, 7.200])
  vol_ref  = np.array([0.865, 0.858, 0.841, 0.811, 0.560, 0.557, 0.546, 0.540])
  str_ref  = np.array([11.43, 12.39, 14.32, 18.67, 44.96, 46.65, 52.13, 55.86])
  vel_ref  = np.array([0.813, 0.867, 0.986, 1.226, 2.907, 2.971, 3.178, 3.311])
  vol_rsd  = np.array([0.006, 0.020, 0.010, 0.005, 0.003, 0.004, 0.004, 0.004])
  str_rsd  = np.array([0.526, 1.708, 0.937, 0.534, 0.280, 0.374, 0.424, 0.437])
  vel_rsd  = np.array([0.037, 0.119, 0.064, 0.035, 0.018, 0.024, 0.026, 0.026])

## Read Data and Calculate Volume ##
  for n in range(ndir):
    vol = np.array([]); vel = np.array([]); str = np.array([]); stp = np.array([])
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
    avol = np.append(avol, np.average(vol[len(stp)-nave-1:len(stp)]/bvol))
    avel = np.append(avel, np.average(vel[len(stp)-nave-1:len(stp)]))
    astr = np.append(astr, np.average(str[len(stp)-nave-1:len(stp)]))
    vol_sd = np.append(vol_sd, np.std(vol[len(stp)-nave-1:len(stp)]/bvol))
    vel_sd = np.append(vel_sd, np.std(vel[len(stp)-nave-1:len(stp)]))
    str_sd = np.append(str_sd, np.std(str[len(stp)-nave-1:len(stp)]))
    
## Output File ## 
  if ldmp:
    fpd = open('avehug.dat', 'w')
    fpd.write('# Shock velocity (mm/us), Average specific volume, Average particle velocity (mm/us), Average Pressure (GPa)\n')
    for i in range (len(vshk)): fpd.write('%s ' %(vshk[i]))
    fpd.write('\n')
    for i in range (len(avol)): fpd.write('%s ' %(avol[i]))
    fpd.write('\n')
    for i in range (len(avel)): fpd.write('%s ' %(avel[i]))
    fpd.write('\n')
    for i in range (len(astr)): fpd.write('%s ' %(astr[i]))
    fpd.write('\n')

## Plot Figure ##
  plt.title('Hugoniot Pressure vs. Specific Volume')
  plt.errorbar(vol_ref, str_ref, xerr=vol_rsd, yerr=str_rsd, capsize=4, marker='o', linestyle='None')  
  #coef = np.polyfit(vol_ref, str_ref, 1)
  #appr = np.poly1d(coef)(vol_ref)
  #plt.plot(vol_ref, appr,  color = 'black', linestyle=':')
  plt.errorbar(avol, astr, xerr=vol_sd, yerr=str_sd, capsize=4, marker='^', linestyle='None')  
  #coef = np.polyfit(x, y, 1)
  #appr = np.poly1d(coef)(x)
  #plt.plot(x, appr,  color = 'black', linestyle=':')
  plt.xlabel('Specific volume')
  plt.ylabel('Pressure (GPa)')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()
  
  plt.title('Shock Velocity vs. Particle Velocity')
  plt.errorbar(vel_ref, vshk_ref, xerr=vel_rsd, capsize=4, marker='o', linestyle='None')  
  plt.errorbar(avel, vshk, xerr=vel_sd, capsize=4, marker='^', linestyle='None', markerfacecolor='None')  
  plt.xlabel('Particle velocity (mm/$\mathrm{\mu}$s)')
  plt.ylabel('Shock velocity (mm/$\mathrm{\mu}$s)')
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()

plot_avHUG()
