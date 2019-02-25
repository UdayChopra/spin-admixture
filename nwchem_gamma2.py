#!/usr/bin/python3

'''
This python module calculates the spin-admixture parameter [citation]. It has
been interfaced with NWChem but can be extended to other Quantum Chemistry codes
with similar inputs. It depends on the 'nwchem_movecs_parser.py' for the wavefunctions
and 'nwchem_zora_parser.py' for the spin-orbit matrix elements.

Input: NWChem wavefunctions (.movecs) and spin-orbit coupling matrix (.zora_so) files.  
Output: Spin-admixture parameter (gamma2). 

Command: % nwchem_gamma2.py movecsfilename zorafilename

Uday Chopra 2018
'''

import argparse
import os
import numpy as np
import nwchem_movecs_parser as nwcparse
import nwchem_zora_parser as zora
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='Uses spin-polarized wavefunctions of scalar relativistic DFT and zora spin-orbit matrix from a fully relativistic DFT calculation to calculate gamma^2')
parser.add_argument("movecsfile", metavar="NWChem .movecs file")
parser.add_argument("zorafile", metavar="NWChem .zora_so file")
parser.add_argument("-q", "--quiet", dest="q", action="store_true",\
	help="Only print g^2 on stdout")
parser.add_argument("-v", "--verbose", dest="v", action="store_true",\
	help="Print cummulative gamma2 for all the orbitals ")

args = parser.parse_args()
outfilename = os.path.splitext(args.movecsfile)[0]

m      = nwcparse.read(args.movecsfile) # .movecs file object

n = 0 # The number of orbitals below HOMO which is to be perturbed. Zero (HOMO) by default for charged pi-conjugated molecules. non-zero in case of d-orbital Metal complexes which usulally lie below HOMO
ha     = m.homo_index[0]- n # index of the perturbed orbital

# alpha and beta spin density-matrix
m.psia = np.matrix(m.psia)      
m.psib = np.matrix(m.psib)      

zora_sf_up, zora_sf_dn, zora_sf_sc_up, zora_sf_sc_dn, zora_so_x, zora_so_y, zora_so_z, zora_so_sc_x, zora_so_sc_y, zora_so_sc_z = zora.read_nwchem_zora(args.zorafile)

nstates = len(m.evalsa)

H_uu = np.array(np.zeros(nstates,dtype=np.complex128))
H_ud = np.array(np.zeros(nstates,dtype=np.complex128))

# Transforming basis of H_soc to spin-polarized basis 
zora_so_x = np.matmul(m.psib.getH(),np.matmul(zora_so_x,m.psia))
zora_so_y = np.matmul(m.psib.getH(),np.matmul(zora_so_y,m.psia))
zora_so_z = np.matmul(m.psia.getH(),np.matmul(zora_so_z,m.psia))

# Deleting the ha (alpha spin HOMO intex) row-vector from the alpha and beta eigen-value differences (perturbation theory doesn't loop over alpha-HOMO)
dea = np.delete((m.evalsa[:] - m.evalsa[ha]),ha)
deb = np.delete((m.evalsb[:] - m.evalsa[ha]),ha)
olx = np.delete(np.asarray(zora_so_x[:,ha]),ha)   # SOC matrix elements <psi_k| H |psi_homo>
oly = np.delete(np.asarray(zora_so_y[:,ha]),ha)
olz = np.delete(np.asarray(zora_so_z[:,ha]),ha)

# Coefficients of the perturbed wavefunction
H_uu = np.float64(0.5) * olz / dea
H_ud = np.float64(0.5) * (olx + 1j*oly) / deb

# up-up and up-down spin-mixing for each orbital
gamma_uu = H_uu.conj()*H_uu     
gamma_ud = H_ud.conj()*H_ud     
gamma2   = gamma_uu + gamma_ud

# Verbose
if (args.v):
    with open("cummulative_"+outfilename+".dat", 'w+') as f:
        f.write(("#Index\t\tgamma2\t\tgamma2_up\t\tgamma2_dn\t\tCumm gamma2\t\tCumm_gamma2_up\t\tCumm_gamma2_dn\n"))
        for i in range(len(H_uu)):
            f.write(("%i\t%16.12e\t%16.12e\t%16.12e\t%16.12e\t%16.12e\t%16.12e\n") %(i+1, gamma2[i].item().real, gamma_uu[i].item().real, gamma_ud[i].item().real, np.sum(gamma2[:i+1]).item().real, np.sum(gamma_uu[:i+1]).item().real, np.sum(gamma_ud[:i+1]).item().real))

gamma_uu = np.sum(gamma_uu)
gamma_ud = np.sum(gamma_ud)

gamma2 = gamma_uu + gamma_ud

if (args.q):
    print(("%16.12e") % gamma2.item().real )
else:
    print("Finished processing .movecs and .zora_so files") 
    fname = outfilename + "_gamma2.dat"
    with open(fname, 'w+') as f:
        f.write(("#Total gamma2\t\tgamma2_up\t\tgamma2_dn\n%16.12e\t%16.12e\t%16.12e\n") %(gamma2.item().real, gamma_uu.item().real, gamma_ud.item().real))

print(("Calculation time : %5.3f minutes") %((time.time()-start_time)/60))

