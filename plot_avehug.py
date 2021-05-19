import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_avehug():
## Number ant Path of the Files ##
  ndir = 5;   dirnam = [[]]*ndir
  dirnam[0] = 'MSST/avehugR.dat'
  dirnam[1] = 'MSST/avehugE.dat'
  dirnam[2] = 'MSST/avehugEF.dat'
  dirnam[3] = 'MSST/avehugEP.dat'
  dirnam[4] = 'MSST/avehugEFP.dat'

## Computational Conditions ##  
  lplot = 2           # 1: P-V  2: Us-Up
  lerr = False         # draw error bar 
  lapp = True         # draw approximate line
  v0   = -23273.904371 # Initial volume (< 0: volume in bohr^3,  >=0: referance step)
  cl = ['k', 'r', 'b', 'g', 'y']
  mk = ['o', '^', 'v', 's', 'x']
  lg = ['DFT', 'E', 'EF', 'EP', 'EFP']
  fc = ['k', 'None', 'None', 'None', 'None']

####################################################################################################

## Read Data and Calculate Volume ##
  for n in range(ndir):
    elpl = 1
    fp = open('%s' %(dirnam[n]), 'r'); fp.readline()
    Us = np.array(fp.readline().split(), dtype='float')
    V  = np.array(fp.readline().split(), dtype='float')
    Up = np.array(fp.readline().split(), dtype='float')
    P  = np.array(fp.readline().split(), dtype='float')
    sdV  = np.array(fp.readline().split(), dtype='float')
    sdUp = np.array(fp.readline().split(), dtype='float')
    sdP  = np.array(fp.readline().split(), dtype='float')
    for i in range(len(Us)):
      if V[i] < 0.7: etp = i; break 
    if lplot == 1:      
      if lerr:
        plt.errorbar(V[:etp], P[:etp], xerr=sdV[:etp], yerr=sdP[:etp], capsize=4, color=cl[n], marker=mk[n], linestyle='None', label=lg[n])  
        plt.errorbar(V[etp:], P[etp:], xerr=sdV[etp:], yerr=sdP[etp:], capsize=4, color=cl[n], marker=mk[n], linestyle='None')  
      else:
        plt.scatter(V[:etp], P[:etp], color=cl[n], marker=mk[n], linestyle='None', label=lg[n])  
        plt.scatter(V[etp:], P[etp:], color=cl[n], marker=mk[n], linestyle='None')  
      if lapp:
        coef = np.polyfit(V[:etp], P[:etp], 1)
        appr = np.poly1d(coef)(V[:etp])
        plt.plot(V[:etp], appr,  color = cl[n], linestyle='--')
        coef = np.polyfit(V[etp:], P[etp:], 1)
        appr = np.poly1d(coef)(V[etp:])
        plt.plot(V[etp:], appr,  color = cl[n], linestyle='--')
    if lplot == 2:      
      if lerr:
        plt.errorbar(Up[:etp], Us[:etp], yerr=sdUp[:etp], capsize=4, color=cl[n], marker=mk[n], linestyle='None', label=lg[n])  
        plt.errorbar(Up[etp:], Us[etp:], yerr=sdUp[etp:], capsize=4, color=cl[n], marker=mk[n], linestyle='None')  
      else:
        plt.scatter(Up[:etp], Us[:etp], color=cl[n], marker=mk[n], linestyle='None', label=lg[n])  
        plt.scatter(Up[etp:], Us[etp:], color=cl[n], marker=mk[n], linestyle='None')  
      if lapp:
        coef = np.polyfit(Up[:etp], Us[:etp], 1)
        appr = np.poly1d(coef)(Up[:etp])
        plt.plot(Up[:etp], appr,  color = cl[n], linestyle='--')
        coef = np.polyfit(Up[etp:], Us[etp:], 1)
        appr = np.poly1d(coef)(Up[etp:])
        plt.plot(Up[etp:], appr,  color = cl[n], linestyle='--')
    
  if lplot == 1:  
    plt.xlabel('Specific volume'); plt.ylabel('Pressure (GPa)')
    plt.title('Hugoniot Pressure vs. Specific Volume')
  if lplot == 2:  
    plt.xlabel('Particle velocity (mm/$\mathrm{\mu}$s)'); plt.ylabel('Shock velocity (mm/$\mathrm{\mu}$s)')
    plt.title('Hugoniot Pressure vs. Specific Volume')

  plt.legend()
  plt.subplots_adjust(left=0.165, right=0.94, bottom=0.12, top=0.92)
  plt.show()

plot_avehug()
