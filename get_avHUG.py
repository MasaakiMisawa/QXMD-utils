import numpy as np

def get_avHUG():
## Number of Directory and Directory Name ##
  ndir = 16;   dirnam = [[]]*ndir
  for n in range(ndir): dirnam[n] = './MSST/%d' %(6000+100*n) 

## Computational Conditions ##  
  nini = 500           # Initial step 
  nend = 1e9           # Final step
  nave = 1500          # Average range [step]
  cskp = 1             # Skip step in md_cel.d
  v0   = -23273.904371 # Initial volume (< 0: volume in bohr^3,  >=0: referance step)

####################################################################################################
# Files to be read: (1) md_cel.d (2) md_hug.d
####################################################################################################

## Constant and Array ##
  lini = 1; bt = 0; bvol = -v0;  
  vol_sd = np.array([]); vel_sd = np.array([]); str_sd = np.array([])
  avol = np.array([]); avel = np.array([]); astr = np.array([])
  vshk = np.array([6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5])

## Read Data and Calculate Volume ##
  for n in range(ndir):
    vol = np.array([]); vel = np.array([]); str = np.array([]); stp = np.array([])
    fpc = open('%s/md_cel.d' %(dirnam[n]), 'r')
    for i in range(2): fpc.readline().split() 
    fph = open('%s/md_hug.d' %(dirnam[n]), 'r')
    for i in range(4): fph.readline().split() 
    while True:
      dth = np.array(fph.readline().split(), dtype='float')
      if len(dth) == 0: break
      if dth[0]%cskp == 0: dtc = np.array(fpc.readline().split(), dtype='float')   
      if len(dtc) == 0: break
      if int(dth[0]) >= nini and int(dth[0]) <= nend: 
        cel = np.array(dtc[1:10]).reshape(3,3).T
        vol0 = np.dot(cel[:,0],np.cross(cel[:,1],cel[:,2])); vol = np.append(vol, vol0)
        stp = np.append(stp, dth[0]); 
        vel = np.append(vel, dth[1])
        str = np.append(str, dth[2])
        if v0 == int(dtc[0]): 
          bvol = vol0
    fpc.close()
    fph.close()
    if len(stp) - nave + 1 <= 0:
      print('Error: too small nave'); return 1
    vel *= 1.e-3  # m/s -> km/s
    avol = np.append(avol, np.average(vol[len(stp)-nave-1:len(stp)]/bvol))
    avel = np.append(avel, np.average(vel[len(stp)-nave-1:len(stp)]))
    astr = np.append(astr, np.average(str[len(stp)-nave-1:len(stp)]))
    vol_sd = np.append(vol_sd, np.std(vol[len(stp)-nave-1:len(stp)]/bvol))
    vel_sd = np.append(vel_sd, np.std(vel[len(stp)-nave-1:len(stp)]))
    str_sd = np.append(str_sd, np.std(str[len(stp)-nave-1:len(stp)]))
    
## Output File ## 
  fpd = open('avehug.dat', 'w')
  fpd.write('# <Us> (mm/us), <v/v0>, <Up> (mm/us), <P> (GPa), SD(Us), SD(V/V0), SD(Up), SD(P) \n')
  for i in range (len(vshk)): fpd.write('%s ' %(vshk[i]))
  fpd.write('\n')
  for i in range (len(avol)): fpd.write('%s ' %(avol[i]))
  fpd.write('\n')
  for i in range (len(avel)): fpd.write('%s ' %(avel[i]))
  fpd.write('\n')
  for i in range (len(astr)): fpd.write('%s ' %(astr[i]))
  fpd.write('\n')
  for i in range (len(vol_sd)): fpd.write('%s ' %(vol_sd[i]))
  fpd.write('\n')
  for i in range (len(vel_sd)): fpd.write('%s ' %(vel_sd[i]))
  fpd.write('\n')
  for i in range (len(str_sd)): fpd.write('%s ' %(str_sd[i]))
  fpd.write('\n')

get_avHUG()
