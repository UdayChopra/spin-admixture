# spin-admixture

This python module calculates the spin-admixture parameter. It has
been interfaced with NWChem but can be extended to other Quantum Chemistry codes
with similar inputs. It depends on the 'nwchem_movecs_parser.py' for the wavefunctions
and 'nwchem_zora_parser.py' for the spin-orbit matrix elements.


# Citation

U. Chopra, Shambhawi, Sergei A. Egorov, Jairo Sinova, Erik R. McNellis "Accurate and general formalism for spin-mixing parameter calculations" Phys. Rev. B. (2019)

# Usage

Input: NWChem wavefunctions (.movecs) and spin-orbit coupling matrix (.zora_so) files.  
Output: Spin-admixture parameter (gamma2). 

Command: % nwchem_gamma2.py movecsfilename zorafilename
