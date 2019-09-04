#!/usr/bin/env python3

'''
This parser reads the NWChem ".zora_so" file in unformatted FORTRAN binary format and returns all the relevant data as 2D numpy matrices.

Uday Chopra 2018
'''

from scipy.io import FortranFile
import numpy as np


def read_ints(f,dt):
	r = np.fromfile(f, dtype=np.int32, count=1)[0]
	r = int(r / np.dtype(dt).itemsize)
	ints = np.fromfile(f, dtype=dt, count=r)
	f.seek(4,1)
	return ints

def read_nwchem_zora(fname):
    with open(fname) as ff:
        l = np.fromfile(ff, dtype=np.int32, count=1)[0]
    itype = np.int32
    if ( l == 8 ): itype = np.int64
    f = FortranFile(fname, 'r')
    nsets = f.read_record(dtype=itype)[0]   # 2 for UHF 1 for RHF
    nbf   = f.read_record(dtype=itype)[0]   # number of basis functions
    mult  = f.read_record(dtype=itype)[0]   # spin-multiplicity

# Variable names of ten (nbf x nbf) matrices are self-explanatory. Only spin-orbit matrices need to be imaginary. 
    
    zora=[]
    zora_sf_up=[]           # spin-orbit free; scalar relativistic parts                       
    zora_sf_dn=[]           # spin-orbit free; scalar relativistic parts           
    zora_sf_sc_up=[]        # spin-orbit free; scaled         
    zora_sf_sc_dn=[]        # spin-orbit free; scaled 
    zora_so_x=[]            # gamma2 uses only these three SOC matrices  
    zora_so_y=[]            # gamma2 uses only these three SOC matrices
    zora_so_z=[]            # gamma2 uses only these three SOC matrices
    zora_so_sc_x=[]         # scaled zora SOC matrix-elements         
    zora_so_sc_y=[]         # scaled zora SOC matrix-elements
    zora_so_sc_z=[]         # scaled zora SOC matrix-elements

    for j in range(nbf):
        zora_sf_up=np.append(zora_sf_up,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_sf_dn=np.append(zora_sf_dn,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_sf_sc_up=np.append(zora_sf_sc_up,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_sf_sc_dn=np.append(zora_sf_sc_dn,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_so_z=np.append(zora_so_z,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_so_y=np.append(zora_so_y,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_so_x=np.append(zora_so_x,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_so_sc_z=np.append(zora_so_sc_z,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_so_sc_y=np.append(zora_so_sc_y,f.read_record(dtype=np.float64))

    for j in range(nbf):
        zora_so_sc_x=np.append(zora_so_sc_x,f.read_record(dtype=np.float64))

    zora_so_x=1j*np.matrix(np.reshape(zora_so_x,(nbf,nbf)))
    zora_so_y=1j*np.matrix(np.reshape(zora_so_y,(nbf,nbf)))
    zora_so_z=1j*np.matrix(np.reshape(zora_so_z,(nbf,nbf))) 
    
    zora = {"sf_up":zora_sf_up,
            "sf_dn":zora_sf_dn,
            "sf_up_sc":zora_sf_sc_up,
            "sf_dn_sc":zora_sf_sc_dn,
            "so_z":zora_so_z,
            "so_y":zora_so_y,
            "so_x":zora_so_x,
            "so_z_sc":zora_so_sc_z,
            "so_y_sc":zora_so_sc_y,
            "so_x_sc":zora_so_sc_x,
            }

    return zora
