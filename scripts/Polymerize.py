import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import re
from mendeleev import element
import sys

def read_input(in_fname):
    f = open(in_fname, "r")
    fl = f.readlines()
    f.close()

    VAR = {}
    for line in fl:
        line = line.split()
        VAR[line[0]] = line[2]

    VAR['LINKER_NDXS'] = np.array(VAR['LINKER_NDXS'].split(','), dtype='int')
    VAR['N_LINKERS'] = np.array(VAR['N_LINKERS'].split(','), dtype='int')
    VAR['CAP_NDX'] = int(VAR['CAP_NDX'])
    VAR['LINKER_DIR'] = str(VAR['LINKER_DIR'])
    VAR['CAP_DIR'] = str(VAR['CAP_DIR'])
    VAR['OUTNAME'] = str(VAR['OUTNAME'])

    return VAR

def read_ATOM(lines):
    read = False
    block = []

    for line in lines:
        if "@<TRIPOS>BOND" in line:
            break
        elif read:
            block.append(line.split())
        elif "@<TRIPOS>ATOM" in line:
            read = True
    block = np.array(block, dtype='str')

    names = block[:,1]
    xyz = block[:,2:5].astype('float')
    atypes = block[:,5]
    charges = block[:,8].astype('float')

    return names, xyz, atypes, charges

def read_BOND(lines):
    read = False

    block = []
    for line in lines:
        if "@<TRIPOS>SUBSTRUCTURE" in line:
            break
        elif read:
            block.append(line.split())
        elif "@<TRIPOS>BOND" in line:
            read = True
    block = np.array(block, dtype='str')

    bonds = block[:,1:3].astype('int')
    btypes = block[:,3]
    return bonds, btypes

def read_RESIDUECONNECT(lines):
    read = False

    block = []
    for line in lines:
        if read:
            block.append(line.split())
        elif "@<TRIPOS>RESIDUECONNECT" in line:
            read = True
    block = np.array(block, dtype='str')

    return block[0]

def rot_mat(u, t):
    #Rotation matrix around vector u for an angle t
    ct = np.cos(t)
    st = np.sin(t)
    x = u[0]
    y = u[1]
    z = u[2]
    rot = np.array([[ct + x**2*(1-ct), x*y*(1-ct)-z*st, x*z*(1-ct)+y*st],     [x*y*(1-ct)+z*st, ct+y**2*(1-ct), y*z*(1-ct)-x*st],    [x*z*(1-ct)-y*st, y*z*(1-ct)+x*st, ct+z**2*(1-ct)]])
    return rot

def weighted_pca(frag, weight=False, noh=True):
    if weight:
        good_xyz = []
        elements = []
        for name in frag.names:
            elements.append(re.split('(\d+)',name)[0])
        for e, el in enumerate(elements):
            for j in range(element(el).atomic_number):
                good_xyz.append(frag.xyz[e])
        good_xyz = np.array(good_xyz)

    elif noh:
        heavy = heavy_ndxs(frag)
        good_xyz = frag.xyz[heavy]
        if np.sum(heavy) == 2:
            good_xyz = np.append(heavy_xyz, [0,0,0])
    else:
        good_xyz = frag.xyz*1

    pca = PCA(n_components=3)
    pca.fit(good_xyz)
    pca1 = pca.components_[0]
    return pca1

def clean_coordinates(frag):
    frag.xyz = frag.xyz - np.mean(frag.xyz, axis=0)

    #ROTATION ONTO PLANE XY
    pca1 = weighted_pca(frag)
    proj_xy = pca1
    proj_xy[2] = 0
    cos = np.dot(pca1, proj_xy)/np.linalg.norm(proj_xy)
    if cos - 1 < 10**-5 and cos > 1:
        cos = 1
    theta = np.arccos(cos)
    ax_rot = np.cross(pca1, proj_xy)
    R = rot_mat(ax_rot, theta)
    frag.xyz = np.dot(R, frag.xyz.T).T

    #ALIGNMENT TO X
    pca1 = weighted_pca(frag)
    proj_x = pca1*1 #The one is to avoid overwriting pca1 in the next line
    proj_x[1:] = [0,0]
    cos = np.dot(pca1, proj_x)/(np.linalg.norm(proj_x))
    if cos - 1 < 10**-5 and cos > 1:
        cos = 1
    theta = np.arccos(cos)
    if (pca1[0] > 0 and pca1[1] > 0) or (pca1[0] < 0 and pca1[1] < 0):
        theta = -theta
    ax_rot = [0,0,1]
    R = rot_mat(ax_rot, theta)
    frag.xyz = np.dot(R, frag.xyz.T).T

    #HEAD-TAIL ORIENTATION
    if frag.xyz[frag.tail,0] > 0:
        R = rot_mat([0,0,1], np.pi)
        frag.xyz = np.dot(R, frag.xyz.T).T

    return frag.xyz

def heavy_ndxs(frag):
    heavy = []
    for name in frag.names:
        if "H" != name[0] and "h" != name[0]:
            heavy.append(True)
        else:
            heavy.append(False)
    return heavy

def find_bonded_ndxs(atid, bonds):
    bonded_ndxs = []
    for bond in bonds:
        if atid == bond[0]:
            bonded_ndxs.append(bond[1] - 1)
        elif atid == bond[1]:
            bonded_ndxs.append(bond[0] - 1)
    return bonded_ndxs

def place_summit(central, bonded):
    bond_vecs = bonded - central
    for i in range(len(bond_vecs)):
        bond_vecs[i] = bond_vecs[i]/np.linalg.norm(bond_vecs[i])
    summit_dir = np.mean(bond_vecs, axis=0)
    summit_dir = -summit_dir/np.linalg.norm(summit_dir)
    return summit_dir

def update_names(frag1, frag2):
    all_names = np.append(frag1.names, frag2.names)
    elements = []
    counts = {}
    new_names = []
    for name in all_names:
        elements.append(re.split('(\d+)',name)[0])
    for el in np.unique(elements):
        counts[el] = 0
    for el in elements:
        counts[el] += 1
        new_names.append("{}{}".format(el, counts[el]))
    new_names = np.array(new_names)
    return new_names

def search_best_orientation(main_xyz, frag, delta, head):
    heavy = heavy_ndxs(frag)
    anchored_xyz = frag.xyz - frag.xyz[frag.tail]
    delta = delta - np.mean(anchored_xyz, axis=0)
    pca = PCA(n_components=3)
    pca.fit(anchored_xyz[heavy])
    pca1 = pca.components_[0]

    psis = np.linspace(0, 2*np.pi, 60)
    best_dist = 0.0
    best_psi = 0.0

    for psi in psis:
        R = rot_mat(pca1, psi)
        test_xyz = np.dot(R, anchored_xyz.T).T
        test_xyz = test_xyz + delta
        dists = cdist(main_xyz, test_xyz)
        dists[head,frag.tail] = 10.0 #Overrides the head-tail distance which would otherwise be always the smallest
        test_dist = np.min(dists)
        if test_dist > best_dist:
            best_dist = test_dist*1
            best_psi = psi*1
    best_R = rot_mat(pca1, best_psi)
    best_xyz = np.dot(best_R, anchored_xyz.T).T
    best_xyz = best_xyz + delta
    return best_xyz

class Fragment:
    def __init__(self, mol2_fname):
        f = open(mol2_fname, "r")
        fl = f.readlines()
        f.close()

        self.names, self.xyz, self.atypes, self.charges = read_ATOM(fl)
        self.bonds, self.btypes = read_BOND(fl)
        terminal = read_RESIDUECONNECT(fl)
        self.tail = np.where(self.names==terminal[1])[0][0]
        if terminal[2] == "0":
            self.head = -1
            self.xyz = clean_coordinates(self)
            self.summit_dir = []
        else:
            self.head = np.where(self.names==terminal[2])[0][0]
            self.xyz = clean_coordinates(self)
            head_bonded_ndxs = find_bonded_ndxs(self.head+1, self.bonds)
            self.summit_dir = place_summit(self.xyz[self.head], self.xyz[head_bonded_ndxs])

    def polymerize(self, frag):
        N_at_prev = len(self.names)

        self.names = update_names(self, frag)

        delta_xyz = self.xyz[self.head] + 1.54*self.summit_dir - frag.xyz[frag.tail]
        frag.xyz = search_best_orientation(self.xyz, frag, delta_xyz, self.head)
        self.xyz = np.vstack((self.xyz, frag.xyz))

        self.atypes = np.append(self.atypes, frag.atypes)
        self.charges = np.append(self.charges, frag.charges)

        self.bonds = np.vstack((self.bonds, frag.bonds+N_at_prev))
        self.bonds = np.vstack((self.bonds, [[self.head+1, frag.tail+N_at_prev+1]]))
        self.btypes = np.append(self.btypes, frag.btypes)
        self.btypes = np.append(self.btypes, "1")

        if frag.head == -1:
            self.head = -1
            self.xyz = clean_coordinates(self)
            self.summit_dir = []
        else:
            self.head = N_at_prev + frag.head
            self.xyz = clean_coordinates(self)
            head_bonded_ndxs = find_bonded_ndxs(self.head+1, self.bonds)
            self.summit_dir = place_summit(self.xyz[self.head], self.xyz[head_bonded_ndxs])

    def write_mol2(self, fname):
        tail_name = self.names[self.tail]
        if self.head == -1:
            head_name = 0
            head_res = 0
        else:
            head_name = self.names[self.head]
            head_red = 0
        N_at, N_bonds = len(self.names), len(self.bonds)

        f = open(fname, "w")
        f.write("@<TRIPOS>MOLECULE\n")
        f.write("THI\n")
        f.write("{}\t\t{}\t\t{}\n".format(N_at, N_bonds, 1))
        f.write("SMALL\nUSER_CHARGES\n\n")

        f.write("@<TRIPOS>ATOM\n")
        atids = np.linspace(1, N_at, N_at, dtype='int')

        for atid, name, coord, atype, charge in zip(atids, self.names, self.xyz, self.atypes, self.charges):
            f.write("{:<3d} {:<5} {:>10.5f} {:>10.5f} {:>10.5f} {:<3} {:>2d} {:>5} {:>7.4f} {:>7.4f} ****\n".format(atid, name, *coord, atype, 1, "THI", charge, 0))

        bids = np.linspace(1, N_bonds, N_bonds, dtype='int')
        f.write("@<TRIPOS>BOND\n")
        for bid, bond, btype in zip(bids, self.bonds, self.btypes):
            f.write("{:>3d} {:>5d} {:>5d} {:<5}\n".format(bid, *bond, btype))

        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write("{:>4d} {:<4} {:>4d} {:<4} {:>4d} {:<4} {:<4}\n".format(1, "THI", 1, "****", 0, "****", "****"))

        f.write("@<TRIPOS>HEADTAIL\n")
        f.write("{:<2} {:<2d}\n{:<2} {:<2d}\n".format(tail_name, 1, head_name, head_res))
        f.write("@<TRIPOS>RESIDUECONNECT\n")
        f.write("{:<1d} {:<2} {:<2} 0 0 0 0\n".format(1, tail_name, head_name))
        f.close()

if __name__ == "__main__":
    VAR = read_input(sys.argv[1])

    thiol = Fragment("{}/FL{}.mol2".format(VAR['LINKER_DIR'], VAR['LINKER_NDXS'][0]))
    VAR['N_LINKERS'][0] -= 1

    for i, link in enumerate(VAR['LINKER_NDXS']):
        for j in range(VAR['N_LINKERS'][i]):
            #The fragment must be initialized iteratively because some attributes change during the polymerization
            next_frag = Fragment("{}/FL{}.mol2".format(VAR['LINKER_DIR'], link))
            thiol.polymerize(next_frag)

    cap = Fragment("{}/FC{}.mol2".format(VAR['CAP_DIR'], VAR['CAP_NDX']))
    thiol.polymerize(cap)
    
    thiol.write_mol2(VAR['OUTNAME'])
