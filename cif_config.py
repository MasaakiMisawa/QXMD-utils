#cif_config.py ver.2021.11.05
import numpy as np; import datetime

def make_config():
## Setup ##

  filename = './crystal.cif'    # filename of CIF to be read
  mul = [2, 2, 2]               # supercell to be created = mul[0]*mul[1]*mul[2] of unit cell
  site_type = [1, 1, 2]         # atom types of each atom site
  type_name = ['Fe','O']        # atom names of each atom type

# Note #################################################################################
# This script generates configuration file for QXMD from CIF.
# An .xyz file is also generated to visualize the geometry in VMD.
# The latest version will be uploaded on "https://github.com/MasaakiMisawa/QXMD-utils".
########################################################################################

## Prepare ##
  la = 0; lb = 0; lc = 0; alp = 0; bet = 0; gam = 0
  rf = open(filename, 'r')
  wf = open('Config.dat', 'w')
  xyz = open('real.xyz', 'w')

## Read CIF ##
  symop = []; lsite = 0; site = []
  while True:
    dat = rf.readline().split()
    if len(dat) == 0: break
    if dat[0] == '_cell_length_a': la = float(dat[1])
    elif dat[0] == '_cell_length_b': lb = float(dat[1])
    elif dat[0] == '_cell_length_c': lc = float(dat[1])
    elif dat[0] == '_cell_angle_alpha': alp = float(dat[1])
    elif dat[0] == '_cell_angle_beta': bet = float(dat[1])
    elif dat[0] == '_cell_angle_gamma': gam = float(dat[1])
    elif dat[0] == '_space_group_symop_operation_xyz':
      while True:
        j = 0; k = 0
        dat = rf.readline().split()
        if len(dat) == 0: break
        if dat[0][0] != '\'': break
        sign = 1; ope = [0, 0, 0, 0.0]
        for i in dat[0]:
          if i == '-': sign = -1
          elif i == 'x': 
            ope[0] = sign; sign = 1
          elif i == 'y': 
            ope[1] = sign; sign = 1
          elif i == 'z':
            ope[2] = sign; sign = 1
          elif i == '/':
            ope[3] = sign*float(dat[0][k-1])/float(dat[0][k+1]); sign = 1
          elif i == ',' or i == '\'':
            if ope != [0, 0, 0, 0.0]:
              symop.append(ope)
              ope = [0, 0, 0, 0.0]; sign = 1
          k = k + 1
        k = 0
    if dat[0][0:11] == '_atom_site_' and lsite == 0:
      j = 0; k = 0
      while True:
        if dat[0] == '_atom_site_fract_x': k = j 
        elif dat[0][0:11] == '_atom_site_': j = j + 1
        elif len(dat) > 2:
          if k == 0: break
          else:
            site.append([float(dat[k]), float(dat[k+1]), float(dat[k+2])])
            lsite = 1
        else: break
        dat = rf.readline().split()
        if len(dat) == 0: break
    if len(symop) > 0 and len(site) > 0: break
  #print(symop)
  #print(site)

## Get cell vectors ##
  pi = np.arccos(-1); rtd = 180/pi;
  print('unit cell: %.3f %.3f %.3f %.3f %.3f %.3f' %(la, lb, lc, alp, bet, gam))
  la = la*mul[0]; lb = lb*mul[1]; lc = lc*mul[2]
  print('supercell: %.3f %.3f %.3f %.3f %.3f %.3f' %(la, lb, lc, alp, bet, gam))
  lx = la                      # inverse relationship in https://docs.lammps.org/Howto_triclinic.html
  xy = lb*np.cos(gam/rtd)
  xz = lc*np.cos(bet/rtd)
  ly = np.sqrt(lb*lb - xy*xy)
  yz = (lb*lc*np.cos(alp/rtd) - xy*xz)/ly
  lz = np.sqrt(lc*lc - xz*xz - yz*yz)
  cel =np.array([lx, 0, 0, xy, ly, 0, xz, yz, lz]); cel = cel.reshape(3,3)
  print('cel vectors t(a b c):')
  print('%10.6f %10.6f %10.6f' %(cel[0][0], cel[0][1], cel[0][2]))
  print('%10.6f %10.6f %10.6f' %(cel[1][0], cel[1][1], cel[1][2]))
  print('%10.6f %10.6f %10.6f' %(cel[2][0], cel[2][1], cel[2][2]))

## Get fractional coordinates ##
  nmax = int(len(symop)*len(site)/3)*mul[0]*mul[1]*mul[2]
  fc = [[] for i in range(nmax)]; tmpf = [0, 0, 0]
  cnt = 0; tlist = []
  for i in range(len(site)):
    for j in range(0, len(symop), 3):
      tmpf[0] = float(site[i][0])*symop[j][0] + float(site[i][1])*symop[j][1] + float(site[i][2])*symop[j][2] + symop[j][3]
      tmpf[1] = float(site[i][0])*symop[j+1][0] + float(site[i][1])*symop[j+1][1] + float(site[i][2])*symop[j+1][2] + symop[j+1][3]
      tmpf[2] = float(site[i][0])*symop[j+2][0] + float(site[i][1])*symop[j+2][1] + float(site[i][2])*symop[j+2][2] + symop[j+2][3] 
      for k in range(3):
        while tmpf[k] >= 1.0 or tmpf[k] < 0.0:
          if tmpf[k] >= 1.0: tmpf[k] = tmpf[k] - 1.0 
          elif tmpf[k] < 0.0: tmpf[k] = tmpf[k] + 1.0 
      lolap = 0
      for k in range(cnt):
        if fc[k] != 0:
          dx = abs(fc[k][0] - tmpf[0]); dy = abs(fc[k][1] - tmpf[1]); dz = abs(fc[k][2] - tmpf[2])
          if dx >= 0.5: dx = 1.0 - dx
          if dy >= 0.5: dy = 1.0 - dy
          if dz >= 0.5: dz = 1.0 - dz
          dr = np.sqrt(dx*dx + dy*dy + dz*dz)
          if dr < 0.5/np.average([la, lb, lc]): lolap = 1
      if lolap == 0:
        for iz in range(mul[0]):
          for iy in range(mul[1]):
            for ix in range(mul[2]):
              fc[cnt] = [tmpf[0]+ix, tmpf[1]+iy, tmpf[2]+iz] 
              tlist.append(site_type[i])
              cnt = cnt + 1
  #print(fc)
  #print(tlist)
  ntot = cnt 
  print('Number of atoms: %d' %ntot)

## Output ##
  dd = datetime.datetime.now()
  wf.write(' %d\n' %ntot); xyz.write(' %d\n' %ntot); xyz.write('Generated by cif_condif.py on %s\n' %dd)
  for i in range(ntot):
    fx = fc[i][0]/mul[0]; fy = fc[i][1]/mul[1]; fz = fc[i][2]/mul[2]
    rx = np.dot(np.array([fx, fy, fz]),cel[:,0]); ry = np.dot(np.array([fx, fy, fz]),cel[:,1]); rz = np.dot(np.array([fx, fy, fz]),cel[:,2])
    wf.write('%2d %.6f %.6f %.6f\n' %(tlist[i], fx, fy, fz))
    xyz.write('%s %.6f %.6f %.6f\n' %(type_name[tlist[i]-1], rx, ry, rz))
  
  rf.close(); wf.close(); xyz.close()
 
make_config()

