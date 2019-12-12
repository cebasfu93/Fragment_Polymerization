import numpy as np
import sys
import re

def read_atnames(fname):
    f = open(fname, "r")
    fl = f.readlines()
    f.close()

    atnames = []
    read = False
    for line in fl:
        if "@<TRIPOS>BOND" in line:
            break
        if read and len(line.split()) > 2:
            atname = line.split()[1]
            atnames.append(atname)
        elif "@<TRIPOS>ATOM" in line:
            read = True
    atnames = np.array(atnames, dtype='str')
    return atnames

def read_RESIDUECONNECT(fname):
    f = open(fname, "r")
    fl = f.readlines()
    f.close()

    read = False
    block = []
    for line in fl:
        if read:
            block.append(line.split())
        elif "@<TRIPOS>RESIDUECONNECT" in line:
            read = True
    block = np.array(block, dtype='str')

    return block[0]

def find_tail(names, resconnect):
    tail_atid = np.where(names==resconnect[1])[0][0] + 1
    return tail_atid

mol2_fname = sys.argv[1]
atnames = read_atnames(mol2_fname)
resconnect = read_RESIDUECONNECT(mol2_fname)
tail = find_tail(atnames, resconnect)
print(tail)
