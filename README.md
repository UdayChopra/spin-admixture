# spin-admixture

This python module calculates the spin-admixture parameter. It has
been interfaced with NWChem but can be extended to other Quantum Chemistry codes
with similar inputs. It depends on the 'nwchem_movecs_parser.py' for the wavefunctions
and 'nwchem_zora_parser.py' for the spin-orbit matrix elements.


# Citation

"Accurate and general formalism for spin-mixing parameter calculations" 
Uday Chopra, Shambhawi, Sergei A. Egorov, Jairo Sinova, Erik R. McNellis Phys. Rev. B 100, 134410 (2019)
https://doi.org/10.1103/PhysRevB.100.134410

"Chemical and Structural Trends in the Spin-Admixture Parameter of Organic Semiconductor Molecules"
Uday Chopra, Sergei A. Egorov, Jairo Sinova, Erik R. McNellis J. Phys. Chem. C, 123, 31, 19112 (2019)
https://doi.org/10.1021/acs.jpcc.9b03409

# Usage

Input: NWChem wavefunctions (.movecs) and spin-orbit coupling matrix (.zora_so) files.  
Output: Spin-admixture parameter (gamma2). 

Command: % nwchem_gamma2.py movecsfilename zorafilename
