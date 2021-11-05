import numpy as np; import math; import time

def make_TRJ():
## Number of Directory and Directory Name ##
  ndir = 1;   dirnam = [[]]*ndir
  dirnam[0] = './DATA/data'

## Computational Conditions ##
  nini = 0     # Initial step number (write trj if step >= nini)
  nskp = 1     # Skip step (write trj if (step - nini)%nskp == 0 )
  nend = 1e9   # Final step number (write trj if step < nend) 
  fmax = 100   # Maximum number of frames to be written
  lout = 0     # 0: Do not output "vx, vy, vz"
               # 1: Output atomic velocities as "vx, vy, vz"
               # 2: Output atomic forces as "vx, vy, vz"
               # 3: Output atomic charge as "vx"
               # 4: Output atomic charge and polarization as "vx" and "vy", respectively
  
####################################################################################################
# Files to be read: 
#   (1) md_log (2) md_box.d (3) md_spc.d (4) qm_ion.d (FPMD) or md_ion.d (MD)
#                *This program also can be used for parallelized classical MD
# Files to be created: config.lammpstrj
####################################################################################################

## Constant, Array, and File Object ##
  pi = np.arccos(-1); rtd = 180/pi; bta = 0.529177210903; lini = 1; cnt = 0; stc = -1 
  ity = np.array([]); sti = -1; dion = 'qm'; ncmd = 1; trj = open('config.lammpstrj', 'w') 
  start = time.time()

## Select MD/FPMD and Open LAMMPSTRJ File ##
  for n in range(ndir):   #Directory loop
    fp = open('%s/md_log' %(dirnam[n]), 'r') 
    for i in range(6): 
      dat = fp.readline().split()  
      if dat[0] == 'QM' and dat[5] == '0': ncmd = int(dat[6]); dion = 'md' 
    fp.close()

## Open Data File and Read Header ##
    fpc = open('%s/md_box.d' %(dirnam[n]), 'r')
    for i in range(2): dtc = fpc.readline().split()
    fpi = np.array([]); fps = np.array([])
    for i in range(ncmd):
      if ncmd == 1: fn = ''
      elif ncmd < 10: fn = '%d' %(i)
      elif ncmd < 100: fn = '%02d' %(i)
      else: fn = '%03d' %(i)
      fpi = np.append(fpi,open('%s/%s_ion.d%s' %(dirnam[n],dion,fn), 'r'))
      fps = np.append(fps,open('%s/md_spc.d%s' %(dirnam[n],fn), 'r'))
      dti = fpi[i].readline().split(); dts = fps[i].readline().split()
      dts = fps[i].readline().split(); ity = np.array(dts[1:], dtype = 'int') 
     
    while True:   #Step loop

## Get MD Step, Number of Atoms, Atomic Scaled Coordinate and Species ##
      ion = np.array([]); spc = np.array([]); ntot = 0; psti = sti
      for i in range(ncmd):
        dti = fpi[i].readline().split()
        if len(dti) != 0:
          if dion == 'qm': 
            sti = int(dti[0]); nty = int(dti[1]); 
            nat = np.array(dti[2:nty+2], dtype='int'); ntot = sum(nat)
          if dion == 'md': 
            sti = int(dti[0]); nat = np.array(dti[1], dtype='int'); 
            nat = np.append(nat,0); ntot += sum(nat)
          fac = float(fpi[i].readline().split()[0])
          if sum(nat) == 0: fpi[i].readline()
          for j in range(math.ceil(sum(nat)/3)):
            ion = np.append(ion, np.array(fpi[i].readline().split(), dtype = 'float'))
          dts = fps[i].readline().split()
          if sum(nat) == 0: fps[i].readline()
          for j in range(math.ceil(sum(nat)/36)):
            spc = np.append(spc, np.array(fps[i].readline().split(), dtype = 'int'))
      if len(dti) == 0 or len(dts) == 0: break
      if sti > nend: break
      index = np.array(list(range(ntot)))
      if sti >= nini and (sti - nini)%nskp == 0: 
        if dion == 'md': index = np.argsort(spc)
        typ = get_atom(ity)

## Get Box Bounds ##
      if stc < sti:
        dtc = np.array(fpc.readline().split(), dtype = 'float') 
        if len(dtc) != 0: 
          stc = int(dtc[0])
          if stc <= sti and len(dtc) != 0:
            ax = bta*dtc[1]
            bx = bta*dtc[2]*np.cos(dtc[6]/rtd)
            by = bta*dtc[2]*np.sin(dtc[6]/rtd)
            cx = bta*dtc[3]*np.cos(dtc[5]/rtd)
            cy = (bta*dtc[2]*bta*dtc[3]*np.cos(dtc[4]/rtd) - bx*cx)/by
            cz = np.sqrt(bta*bta*dtc[3]*dtc[3] - cx*cx - cy*cy)
            xlo = 0 + min(0.0, bx, cx, bx+cx)
            xhi = ax + max(0.0, bx, cx, bx+cx)
            ylo = 0 + min(0.0, cy)
            yhi = by + max(0.0, cy)
            zlo = 0
            zhi = cz

## Write LAMMPSTRJ File ##
      if sti >= nini and (sti - nini)%nskp == 0:
        if(sti != psti):
          nid = 0; cnt += 1
          trj.write('ITEM: TIMESTEP\n%d\n' %sti)
          trj.write('ITEM: NUMBER OF ATOMS\n%d\n' %ntot)
          trj.write('ITEM: BOX BOUNDS xy xz yz xx yy zz\n')
          trj.write('%.5f %.5f %.5f\n' %(xlo, xhi, bx))
          trj.write('%.5f %.5f %.5f\n' %(ylo, yhi, cx))
          trj.write('%.5f %.5f %.5f\n' %(zlo, zhi, cy))
          trj.write('ITEM: ATOMS element xs ys zs')
          #trj.write(' vx vy vz')
          trj.write('\n')
          for j in index:
            nid += 1      
            xyz = np.array(ion[j*3: j*3 + 3])*fac
            trj.write('%s %7.5f %7.5f %7.5f\n' %(typ[int(spc[j])-1], xyz[0], xyz[1], xyz[2]) )
          if cnt%100 == 0 or lini == 1: print('frame %d (step %d) complete' %(cnt,sti)); lini = 0        
          if fmax <= cnt: break

## Close Files ##
    fpc.close()
    for i in range(ncmd): fpi[i].close(); fps[i].close
    if sti > nend or fmax <= cnt: break
  print('Total number of frames : %d' %(cnt))
  print('Elapsed time: %.2f sec' %(time.time() - start))
  trj.close()

def get_atom(ity):
  typ = np.array([])
  pt = ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
  for i in range (len(ity)): typ = np.append(typ, pt[ity[i]])
  return typ
 
make_TRJ()
