import numpy as np
import sys

def makebonds():
    lammpstrj_name = './config.lammpstrj'
    number_of_types = 3;       cutoff_distances = [[0.0 for i in range(number_of_types)] for j in range(number_of_types)]

    cutoff_distances[0][0] = 0.0
    cutoff_distances[0][1] = 0.0
    cutoff_distances[0][2] = 2.4
    cutoff_distances[1][1] = 0.0
    cutoff_distances[1][2] = 2.1
    cutoff_distances[2][2] = 0.0

    write_lammpstrj = 1

    #################################################################################################

    for i in range(number_of_types - 1):
        for j in range(i+1, number_of_types):
            cutoff_distances[j][i] = cutoff_distances[i][j]

    infile = open('%s' %(lammpstrj_name), 'r')
    outfile = open('bonds.dat', 'w')
    if write_lammpstrj: outfile2 = open('coord.lammpstrj', 'w')
    while (1):  #frame loop
        while (1):
            dat = infile.readline().split()
            if len(dat) == 0:
                sys.exit()
            elif len(dat) == 1: 
                continue
            elif dat[1] == 'TIMESTEP': 
                step = int(infile.readline())
            elif dat[1] == 'NUMBER': 
                ntot = int(infile.readline())
            elif dat[1] == 'BOX':
                dat = infile.readline().split()
                xlob = float(dat[0])
                xhib = float(dat[1])
                xy   = float(dat[2])
                dat = infile.readline().split()
                ylob = float(dat[0])
                yhib = float(dat[1])
                xz   = float(dat[2])
                dat = infile.readline().split()
                zlob = float(dat[0])
                zhib = float(dat[1])
                yz   = float(dat[2])
            elif dat[1] == 'ATOMS': 
                break
        lengths_angles = [0, 0, 0, 0, 0, 0]
        lengths_angles[0] = xhib - max(0, xy, xz, xy+xz)
        lengths_angles[1] = yhib - max(0, yz)
        lengths_angles[2] = zhib
        lengths_angles[3] = deg_rad*np.arccos((xy*xz + np.sqrt(lengths_angles[1]**2 - xy**2)*yz)/(lengths_angles[1]*lengths_angles[2]))  
        lengths_angles[4] = deg_rad*np.arccos(xz/lengths_angles[2])
        lengths_angles[5] = deg_rad*np.arccos(xy/lengths_angles[1])
        cel_vectors = get_celvectors(lengths_angles)
        element = list(range(ntot))
        fractional_coords = list(range(ntot))
        typeindex = []; dic = []
        j = 0
        for i in range(ntot):
            dat = infile.readline().split()
            element[i] = dat[0]
            if len(typeindex) == 0: 
                typeindex.append(j)
                dic.append(dat[0])
                j = j + 1
            else:
                flag = 0
                for k in range(len(dic)):
                    if dic[k] == dat[0]:
                        typeindex.append(k)
                        flag = 1
                        break
                if flag == 0:
                    dic.append(dat[0])
                    typeindex.append(j)
                    j = j + 1
            fractional_coords[i] = [float(dat[1]), float(dat[2]), float(dat[3])]
        #print(cutoff_distances)
        #print(element)
        #print(typeindex)
        #print(len(element))
        #print(len(typeindex))
        #print(fractional_coords)
        #print(lengths_angles)
        #print(cel_vectors)
        neighbor_list = get_neighborlist(fractional_coords, cel_vectors, typeindex, cutoff_distances)
        if write_lammpstrj:
            outfile2.write('ITEM: TIMESTEP\n')
            outfile2.write('%d\n' %step)
            outfile2.write('ITEM: NUMBER OF ATOMS\n')
            outfile2.write('%d\n' %ntot)
            outfile2.write('ITEM: BOX BOUNDS xy xz yz\n')
            outfile2.write('%.5f %.5f %.5f\n' %(xlob, xhib, xy))
            outfile2.write('%.5f %.5f %.5f\n' %(ylob, yhib, xz))
            outfile2.write('%.5f %.5f %.5f\n' %(zlob, zhib, yz))
            outfile2.write('ITEM: ATOMS element xs ys zs vx\n')
        for i in range(ntot):
            #print(neighbor_list[i])
            if write_lammpstrj: outfile2.write('%s %.5f %.5f %.5f %d\n' %(element[i], fractional_coords[i][0], fractional_coords[i][1], fractional_coords[i][2], len(neighbor_list[i])))
            for j in range(len(neighbor_list[i])):
                outfile.write('%d ' %neighbor_list[i][j])
            outfile.write('\n')
    infile.close()
    outfile.close()
    if write_lammpstrj: outfile2.close()

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

def get_neighborlist(fractional_coords, cel_vectors, typeindex, cutoff_distances):
	number_of_atoms = len(typeindex)
	fractional_coords = np.array(fractional_coords)
	neighborlist = [[] for i in range(number_of_atoms)]
	cel_vectors = np.array(cel_vectors, dtype='float')
	if len(fractional_coords) != number_of_atoms:
		print('Error'); return -1
	if cel_vectors.shape != (3, 3):
		print('Error'); return -1
	if np.array(cutoff_distances).shape != (max(typeindex)+1,max(typeindex)+1):
		print('Error'); return -1
	for i in range(number_of_atoms-1):
		for j in range(i+1, number_of_atoms):
			if cutoff_distances[typeindex[i]][typeindex[j]] == 0.0: continue
			fractional_distance = fractional_coords[j] - fractional_coords[i]
			for k in range(3):
				if abs(fractional_distance[k]) > 0.5:
					fractional_distance[k] = fractional_distance[k] - np.sign(fractional_distance[k])
			real_distance = np.dot(cel_vectors.T, fractional_distance)
			if np.linalg.norm(real_distance) <= cutoff_distances[typeindex[i]][typeindex[j]]:
				neighborlist[i].append(j); neighborlist[j].append(i)
	return neighborlist

## Global variables ##
deg_rad = 180.0/np.arccos(-1.0)

makebonds()
