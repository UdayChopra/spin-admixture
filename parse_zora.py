#!/usr/bin/env python3

from scipy.io import FortranFile
import numpy as np

def read_nwchem_zora(fname):
    f = FortranFile(fname, 'r')
    
    nsets=f.read_record(dtype=np.int64)[0]
    nbf=f.read_record(dtype=np.int64)[0]   
    mult=f.read_record(dtype=np.int64)[0]    

# Variable names of ten (nbf x nbf) matrices are self-explanatory. Spin-orbit matrices only need to be imaginary. 
    
    zora=[]
    zora_sf_up=[]
    zora_sf_dn=[]
    zora_sf_sc_up=[]
    zora_sf_sc_dn=[]
    zora_so_x=[]
    zora_so_y=[]
    zora_so_z=[]
    zora_so_sc_x=[]
    zora_so_sc_y=[]
    zora_so_sc_z=[]

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

    zora_so_sc_x=1j*np.matrix(np.reshape(zora_so_sc_x,(nbf,nbf)))
    zora_so_sc_y=1j*np.matrix(np.reshape(zora_so_sc_y,(nbf,nbf)))
    zora_so_sc_z=1j*np.matrix(np.reshape(zora_so_sc_z,(nbf,nbf))) 

    return zora_sf_up, zora_sf_dn, zora_sf_sc_up, zora_sf_sc_dn, zora_so_x, zora_so_y, zora_so_z, zora_so_sc_x, zora_so_sc_y, zora_so_sc_z
    #zora=np.reshape(zora,(nbf,nbf,10))

    #zora_sf_up = zora[:,:,0]
    #zora_sf_dn = zora[:,:,1]
    #zora_sf_sc_up = zora[:,:,2]
    #zora_sf_sc_dn = zora[:,:,3]
    #zora_so_x = np.matrix(1j*zora[:,:,4],dtype=np.complex128)
    #zora_so_y = np.matrix(1j*zora[:,:,5],dtype=np.complex128)
    #zora_so_z = np.matrix(1j*zora[:,:,6],dtype=np.complex128)
    #zora_so_sc_x = np.matrix(1j*zora[:,:,7],dtype=np.complex128)
    #zora_so_sc_y = np.matrix(1j*zora[:,:,8],dtype=np.complex128)
    #zora_so_sc_z = np.matrix(1j*zora[:,:,9],dtype=np.complex128)
