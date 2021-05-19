import numpy as np; import math; import time

def make_PDB():
## Number of Directory and Directory Name ##
  ndir = 1;   dirnam = [[]]*ndir
  dirnam[0] = './MSST/data2'

## Computational Conditions ##
  nini = 0     # Initial step number (write pdb if step >= nini)
  nskp = 1     # Skip step (write pdb if (step - nini)%nskp == 0 )
  nend = 1e9   # Final step number (write pdb if step < nend) 
  fmax = 100   # Maximum number of frames to be written
  
####################################################################################################
# Files to be read: 
#   (1) md_log (2) md_cel.d (3) md_spc.d (4) qm_ion.d (FPMD) or md_ion.d (MD)
#                *This program also can be used for parallelized classical MD
# Files to be created: 
#   config.pdb (including cell lengths and angles for each frame)
####################################################################################################

## Constant, Array, and File Object ##
  pi = np.arccos(-1); rtd = 180/pi; bta = 0.529177210903; lini = 1; cnt = 0; stc = -1 
  ity = np.array([]); sti = -1; dion = 'qm'; ncmd = 1; pdb = open('config.pdb', 'w') 
  start = time.time()

## Read Data and Make PDB File ##
  for n in range(ndir):   #Directory loop
    fp = open('%s/md_log' %(dirnam[n]), 'r') 
    for i in range(6): 
      dat = fp.readline().split()  
      if dat[0] == 'QM' and dat[5] == '0': ncmd = int(dat[6]); dion = 'md' 
    fp.close()

## Open File and Read Header ##
    fpc = open('%s/md_cel.d' %(dirnam[n]), 'r')
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
        if dti != []:
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
      if dti == [] or dts == []: break
      if sti > nend: break
      index = np.array(list(range(ntot)))
      if sti >= nini and (sti - nini)%nskp == 0: 
        if dion == 'md': index = np.argsort(spc)
        typ = get_atom(ity)

## Get Cell Vectors ##
      if stc < sti:
        dtc = fpc.readline().split() 
        if dtc != []: stc = int(dtc[0])

## Get Real Coordinate and Cell Shape and Write PDB File ##
      if sti >= nini and (sti - nini)%nskp == 0:
        if stc <= sti and dtc != []:
          cel = np.array(dtc[1:10], dtype='float').reshape(3,3).T
          a = np.linalg.norm(cel[:,0]); b = np.linalg.norm(cel[:,1]); c = np.linalg.norm(cel[:,2])
          al = np.arccos(np.dot(cel[:,1],cel[:,2])/(b*c))*rtd
          be = np.arccos(np.dot(cel[:,0],cel[:,2])/(a*c))*rtd
          ga = np.arccos(np.dot(cel[:,1],cel[:,0])/(b*a))*rtd
          rmat = 1
          if abs(cel[1,0]) > 1.e-5 or abs(cel[2,0]) > 1.e-5:
            ux = np.array([1,0,0]); 
            ay = np.array([cel[0,0],cel[1,0],0]); az = np.array([cel[0,0],0,cel[2,0]]) 
            db = -np.arccos(np.dot(ux,az)/np.linalg.norm(az))*np.sign(cel[2,0]) 
            dg = np.arccos(np.dot(ux,ay)/np.linalg.norm(ay))*np.sign(cel[1,0])
            rdb = np.matrix(([np.cos(db),0,-np.sin(db)], [0,1,0], [np.sin(db),0,np.cos(db)]))
            rdg = np.matrix(([np.cos(dg),np.sin(dg),0], [-np.sin(dg),np.cos(dg),0], [0,0,1]))
            rmat = rdb*rdg
        if(sti != psti):
          nid = 0; cnt += 1
          pdb.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%16d\n' %(a*bta,b*bta,c*bta,al,be,ga,ntot))
          for j in index:
            nid += 1      
            xyz = np.array(ion[j*3: j*3 + 3])*fac
            xyz = bta*(rmat*np.matrix(cel)*np.matrix(xyz).T).T
            pdb.write('HETATM%5d%3s' %(nid,typ[int(spc[j])-1]))
            pdb.write('%24.3f%8.3f%8.3f%24s\n' %(xyz[0,0],xyz[0,1],xyz[0,2],typ[int(spc[j])-1]))
          pdb.write('END\n')
          if cnt%100 == 0 or lini == 1: print('frame %d (step %d) complete' %(cnt,sti)); lini = 0        
          if fmax <= cnt: break

## Close Files ##
    fpc.close()
    for i in range(ncmd): fpi[i].close(); fps[i].close
    if sti > nend or fmax <= cnt: break
  print('Total number of frames in ./config.pdb: %d' %(cnt))
  print('Elapsed time: %.2f sec' %(time.time() - start))
  pdb.close()

def get_atom(ity):
  typ = np.array([])
  pt = ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
  for i in range (len(ity)): typ = np.append(typ, pt[ity[i]])
  return typ
 
make_PDB()
