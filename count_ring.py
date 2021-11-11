import numpy as np
deg_rad = 180.0/np.arccos(-1.0)

def main():
  directory_name = '../data'
  target_step = 0
  skip_step = 1
  cutoff_distance = [[0, 2.1], [2.1, 0]]

  box = read_box(directory_name, target_step, target_step, skip_step)
  ion = read_ion(directory_name, target_step, target_step, skip_step)
  spc = read_spc(directory_name, target_step, target_step, skip_step)
  cel = get_celvectors(box)
  nlist = get_neighborlist(ion, cel, spc, cutoff_distance)
  glist = get_Gringlist(nlist, ion)
  print('Number of G-rings: %d' %len(glist))
  for i in range(len(glist)): 
    for j in range(len(glist[i])): 
      glist[i][j] = glist[i][j] + 1
  for k in range(20):
    for i in range(len(glist)): 
      if len(glist[i]) == k: print(glist[i])
  print('\n')
  klist = get_Kringlist(nlist, ion)
  print('Number of K-rings: %d' %len(klist))
  for i in range(len(klist)): 
    for j in range(len(klist[i])): 
      klist[i][j] = klist[i][j] + 1
  for k in range(20):
    for i in range(len(klist)): 
      if len(klist[i]) == k: print(klist[i])


def read_box(directory_name, initial_step, final_step, skip_step):
  box = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; box_step = []; typelist =[]
  file = open('%s/md_box.d' %directory_name, 'r')
  bohr = 0.529177210903

  data = file.readline()
  data = file.readline()
  while True:
    data = file.readline().split()
    if len(data) == 0: break
    step = int(data[0])
    if step > final_step: break
    elif step >= initial_step and np.mod(step - initial_step,skip_step) == 0:
      for i in range(3): box[i] = float(data[i+1])*bohr
      for i in range(3): box[i+3] = float(data[i+4])
      box_step.append(box)
  
  if initial_step == final_step: return box_step[0]
  else: return box_step


def read_ion(directory_name, initial_step, final_step, skip_step):
  fractional_coords = []; fractional_coords_step = []
  file = open('%s/qm_ion.d' %directory_name, 'r')

  data = file.readline()
  while True:
    data = np.array(file.readline().split(), dtype = 'int')
    if len(data) == 0: break
    step = data[0]
    number_of_types = data[1]
    number_of_atoms = 0
    for i in range(number_of_types):
      for j in range(data[2+i]):
        number_of_atoms = number_of_atoms + 1
    data = file.readline().split()
    scale_factor = float(data[0])
    if step > final_step: break
    elif step >= initial_step and np.mod(step - initial_step,skip_step) == 0:
      atom_index = 0
      for i in range(int(0.9 + number_of_atoms/3.0)):
        data = np.array(file.readline().split(), dtype = 'float')
        for j in range(3):
          fractional_coords.append([data[np.mod(j,3)*3]*scale_factor, data[np.mod(j,3)*3 + 1]*scale_factor, data[np.mod(j,3)*3 + 2]*scale_factor])
          atom_index = atom_index + 1
          if(atom_index == number_of_atoms): break
      fractional_coords_step.append(fractional_coords)
    else:
      for i in range(int(0.9 + number_of_atoms/3.0)): file.readline()

  if initial_step == final_step: return fractional_coords_step[0]
  else: return fractional_coords_step


def get_neighborlist(fractional_coords, cel_vectors, typeindex, cutoff_distances):
  number_of_atoms = len(typeindex)
  neighborlist = [[] for i in range(number_of_atoms)]
  fractional_distance = np.array([0.0, 0.0, 0.0])
  real_distance = np.array([0.0, 0.0, 0.0])
  cel_vectors = np.array(cel_vectors, dtype='float')

  if len(fractional_coords) != number_of_atoms:
    print('Error1'); return -1 
  if cel_vectors.shape != (3, 3): 
    print('Error2'); return -1 
  if np.array(cutoff_distances).shape != (max(typeindex)+1,max(typeindex)+1):
    print('Error3'); return -1 

  for i in range(number_of_atoms):
    for j in range(number_of_atoms):
      if i == j: continue
      for k in range(3):
        fractional_distance[k] = fractional_coords[j][k] - fractional_coords[i][k]
        if fractional_distance[k] > 0.5: fractional_distance[k] = fractional_distance[k] - 1.0
        if fractional_distance[k] < -0.5: fractional_distance[k] = fractional_distance[k] + 1.0
      for k in range(3):
        real_distance[k] = np.dot(fractional_distance,cel_vectors[:,k])
      if np.linalg.norm(real_distance) <= cutoff_distances[typeindex[i]][typeindex[j]]:
        neighborlist[i].append(j)

  return neighborlist


def get_celvectors(lengths_angles):
  if len(lengths_angles) != 6: 
    print('Error'); return -1 

  ax = lengths_angles[0]
  ay = 0.0
  az = 0.0                     
  bx = lengths_angles[1]*np.cos(lengths_angles[5]/deg_rad)
  by = np.sqrt(lengths_angles[1]*lengths_angles[1] - bx*bx)
  bz = 0.0
  cx = lengths_angles[2]*np.cos(lengths_angles[4]/deg_rad)
  cy = (lengths_angles[1]*lengths_angles[2]*np.cos(lengths_angles[3]/deg_rad) - bx*cx)/by
  cz = np.sqrt(lengths_angles[2]*lengths_angles[2] - cx*cx - cy*cy)
  cel_vectors = [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]

  return cel_vectors


def read_spc(directory_name, initial_step, final_step, skip_step):
  typeindex = []; typeindex_step = []
  file = open('%s/md_spc.d' %directory_name, 'r')

  data = file.readline()
  data = file.readline()
  while True:
    data = np.array(file.readline().split(), dtype = 'int')
    if len(data) == 0: break
    step = data[0]
    number_of_atoms = data[1]
    if step > final_step: break
    elif step >= initial_step and np.mod(step - initial_step,skip_step) == 0:
      atom_index = 0
      while True:
        data = np.array(file.readline().split(), dtype = 'int')
        for j in range(len(data)):
          typeindex.append(data[j]-1)
          atom_index = atom_index + 1
        if(atom_index == number_of_atoms): break
      typeindex_step.append(typeindex)
    else:
      atom_index = 0
      while True:
        data = np.array(file.readline().split(), dtype = 'int')
        for j in range(len(data)):
          atom_index = atom_index + 1
        if(atom_index == number_of_atoms): break

  if initial_step == final_step: return typeindex_step[0]
  else: return typeindex_step


def get_Gringlist(neighborlist, fractional_coords):
  Gringlist = []
  number_of_atoms = len(neighborlist)

  for i in range(number_of_atoms):
    if(len(neighborlist[i]) < 2): continue
    for j in range(len(neighborlist[i])):
      book = []
      book.append([i, neighborlist[i][j]]) 
      while True:
        bookmark = []
        number_of_newpages = 0
        bookpages = len(book)
        for k in range(bookpages):
          for l in range(len(neighborlist[ book[k][len(book[k])-1] ])):
            if neighborlist[book[k][len(book[k])-1]][l] == book[k][len(book[k])-2]: continue
            find_new_member = 1
            if neighborlist[book[k][len(book[k])-1]][l] == book[k][1]:
              n = 1; fractional_distance = [0, 0, 0]; total_fractional_distance = [0,0,0]
              while n+1 < len(book[k]):
                for m in range(3):
                  fractional_distance[m] = fractional_coords[book[k][n+1]][m] - fractional_coords[book[k][n]][m]
                  if abs(fractional_distance[m]) > 0.5: fractional_distance[m] = fractional_distance[m] - np.sign(fractional_distance[m])*1.0
                total_fractional_distance = np.array(total_fractional_distance) + np.array(fractional_distance)
                n = n + 1
              for m in range(3):
                fractional_distance[m] = fractional_coords[book[k][1]][m] - fractional_coords[book[k][n]][m]
                if abs(fractional_distance[m]) > 0.5: fractional_distance[m] = fractional_distance[m] - np.sign(fractional_distance[m])*1.0
              total_fractional_distance = np.array(total_fractional_distance) + np.array(fractional_distance)
              if np.linalg.norm(total_fractional_distance) < 1.0e-5: find_new_member = 0
            if find_new_member == 1: 
              number_of_newpages = number_of_newpages + 1
              newpage = []
              for m in range(len(book[k])):
                newpage.append(book[k][m])
              newpage.append(neighborlist[book[k][len(book[k])-1]][l])
              book.append(newpage)
              newindex = len(book)-1
              if neighborlist[book[k][len(book[k])-1]][l] == i:
                n = 0; fractional_distance = [0, 0, 0]; total_fractional_distance = [0, 0, 0]
                while n+1 < len(book[newindex]):
                  for m in range(3):
                    fractional_distance[m] = fractional_coords[book[newindex][n+1]][m] - fractional_coords[book[newindex][n]][m]
                    if abs(fractional_distance[m]) > 0.5: fractional_distance[m] = fractional_distance[m] - np.sign(fractional_distance[m])*1.0
                  total_fractional_distance = np.array(total_fractional_distance) + np.array(fractional_distance)
                  n = n + 1
                if np.linalg.norm(total_fractional_distance) < 1.0e-5:
                  bookmark.append(newindex)
        if len(bookmark) != 0 or len(book) == bookpages: break
        book = book[bookpages:len(book)]
      if len(bookmark) != 0:
        for k in range(len(bookmark)): 
          ring_member = book[bookmark[k]][0:len(book[bookmark[k]])-1]
          find_new_ring = 1
          for l in range(len(Gringlist)): 
            if Gringlist[l][0] != ring_member[0]:
              for n in range(len(ring_member)):
                ring_member = ring_member[1:len(ring_member)] + [ring_member[0]]
                if Gringlist[l] == ring_member: find_new_ring = 0
            ring_member.reverse()
            ring_member =  [ring_member[len(ring_member)-1]] + ring_member[0:len(ring_member)-1]
            if Gringlist[l] == ring_member: find_new_ring = 0
            if Gringlist[l][0] != ring_member[0]:
              for n in range(len(ring_member)):
                ring_member = ring_member[1:len(ring_member)] + [ring_member[0]]
                if Gringlist[l] == ring_member: find_new_ring = 0
            ring_member.reverse()
            ring_member =  [ring_member[len(ring_member)-1]] + ring_member[0:len(ring_member)-1]
          if find_new_ring == 1: Gringlist.append(ring_member)
  return Gringlist


def get_Kringlist(neighborlist, fractional_coords):
  Kringlist = []
  number_of_atoms = len(neighborlist)

  for i in range(number_of_atoms):
    if(len(neighborlist[i]) < 2): continue
    for j1 in range(len(neighborlist[i])-1):
      for j2 in range(j1,len(neighborlist[i])):
        if len(neighborlist[neighborlist[i][j1]]) < 2 or len(neighborlist[neighborlist[i][j2]]) < 2: continue
        if j1 == j2: continue
        book = []
        book.append([i, neighborlist[i][j1]]) 
        while True:
          bookmark = []
          number_of_newpages = 0
          bookpages = len(book)
          for k in range(bookpages):
            for l in range(len(neighborlist[ book[k][len(book[k])-1] ])):
              if neighborlist[book[k][len(book[k])-1]][l] == book[k][len(book[k])-2]: continue
              find_new_member = 1
              if neighborlist[book[k][len(book[k])-1]][l] == book[k][1] or neighborlist[book[k][len(book[k])-1]][l] == i:
                n = 0
                if neighborlist[book[k][len(book[k])-1]][l] == book[k][1]: n = 1
                fractional_distance = [0, 0, 0]; total_fractional_distance = [0,0,0]
                while n+1 < len(book[k]):
                  for m in range(3):
                    fractional_distance[m] = fractional_coords[book[k][n+1]][m] - fractional_coords[book[k][n]][m]
                    if abs(fractional_distance[m]) > 0.5: fractional_distance[m] = fractional_distance[m] - np.sign(fractional_distance[m])*1.0
                  total_fractional_distance = np.array(total_fractional_distance) + np.array(fractional_distance)
                  n = n + 1
                for m in range(3):
                  fractional_distance[m] = fractional_coords[neighborlist[book[k][len(book[k])-1]][l]][m] - fractional_coords[book[k][n]][m]
                  if abs(fractional_distance[m]) > 0.5: fractional_distance[m] = fractional_distance[m] - np.sign(fractional_distance[m])*1.0
                total_fractional_distance = np.array(total_fractional_distance) + np.array(fractional_distance)
                if np.linalg.norm(total_fractional_distance) < 1.0e-5: find_new_member = 0
              if find_new_member == 1: 
                number_of_newpages = number_of_newpages + 1
                newpage = []
                for m in range(len(book[k])):
                  newpage.append(book[k][m])
                newpage.append(neighborlist[book[k][len(book[k])-1]][l])
                book.append(newpage)
                newindex = len(book)-1
                if neighborlist[book[k][len(book[k])-1]][l] == neighborlist[i][j2]:
                  n = 0; fractional_distance = [0, 0, 0]; total_fractional_distance = [0, 0, 0]
                  while n+1 < len(book[newindex]):
                    for m in range(3):
                      fractional_distance[m] = fractional_coords[book[newindex][n+1]][m] - fractional_coords[book[newindex][n]][m]
                      if abs(fractional_distance[m]) > 0.5: fractional_distance[m] = fractional_distance[m] - np.sign(fractional_distance[m])*1.0
                    total_fractional_distance = np.array(total_fractional_distance) + np.array(fractional_distance)
                    n = n + 1
                  for m in range(3):
                    fractional_distance[m] = fractional_coords[i][m] - fractional_coords[book[newindex][n]][m]
                    if abs(fractional_distance[m]) > 0.5: fractional_distance[m] = fractional_distance[m] - np.sign(fractional_distance[m])*1.0
                  total_fractional_distance = np.array(total_fractional_distance) + np.array(fractional_distance)
                  if np.linalg.norm(total_fractional_distance) < 1.0e-5:
                    bookmark.append(newindex)
          if len(bookmark) != 0 or len(book) == bookpages: break
          book = book[bookpages:len(book)]
        if len(bookmark) != 0:
          for k in range(len(bookmark)): 
            ring_member = book[bookmark[k]][0:len(book[bookmark[k]])]
            find_new_ring = 1
            for l in range(len(Kringlist)): 
              if Kringlist[l][0] != ring_member[0]:
                for n in range(len(ring_member)):
                  ring_member = ring_member[1:len(ring_member)] + [ring_member[0]]
                  if Kringlist[l] == ring_member: find_new_ring = 0
              ring_member.reverse()
              ring_member =  [ring_member[len(ring_member)-1]] + ring_member[0:len(ring_member)-1]
              if Kringlist[l] == ring_member: find_new_ring = 0
              if Kringlist[l][0] != ring_member[0]:
                for n in range(len(ring_member)):
                  ring_member = ring_member[1:len(ring_member)] + [ring_member[0]]
                  if Kringlist[l] == ring_member: find_new_ring = 0
              ring_member.reverse()
              ring_member =  [ring_member[len(ring_member)-1]] + ring_member[0:len(ring_member)-1]
            if find_new_ring == 1: Kringlist.append(ring_member)
  return Kringlist

main()
