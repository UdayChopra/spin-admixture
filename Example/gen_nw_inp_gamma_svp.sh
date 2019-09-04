#!/usr/bin/env bash

#Generates an NWChem input for zora_sf and zora_so calculation reading input from an .XYZ file.

for i in $@; do
	charge=${charge:-"1"}
        xc=${xc:-"pbe0"}
	xcc=${xc}
	[[ ${xc} == "pw92" ]] && xcc="slater pw91lda"
	[[ ${xc} == "pbe" ]] && xcc="xpbe96 cpbe96"
        basis=${basis:-"ma-zora-def2-svp"}
        j1=$(basename -s .xyz ${i})
        j=${j1%}_gamma2_${basis}_${xc}
        mkdir -p scratch
        for task in dft sodft ; do
                infile=${j}_${task}
                cat << EOF > ${infile}.nw
echo

scratch_dir ./scratch

title "${infile}"

Charge ${charge} 

basis "ao basis" spherical
  * library ${basis} except C
  C    S
     1238.40169380             5.727408300E-03
      186.29004992             4.213600060E-02
       42.25117635             1.866435633E-01
       11.67655793             4.678258158E-01
        3.59305065             4.421228295E-01
  C    S
        0.40245147             1.000000000E+00
  C    S
        0.13090183            -1.000000000E+00
  C    P
        9.46809706             5.701225280E-02
        2.01035451             3.130476217E-01
        0.54771005             7.605188875E-01
  C    P
        0.15268614             1.000000000E+00
  C    D
        0.80000000             1.000000000E+00
  #C    S
  #     0.04363394             1.000000000E+00
end

set tolguess 5e-8
set int:acc_std 1d-14

#memory total 1800 mb 

dft
  direct
  xc ${xcc}
  decomp
  grid xfine nodisk     # Extra fine integration grid not stored on disk
  tolerances tight tol_rho 1e-13 accCoul 14     # Coulomb screening tolerances etc.
  convergence energy 5e-9 nolevelshifting       # SCF total energy convergence criterion
  iterations 200        # Maximum SCF iterations
  mult 2
end

relativistic
  zora on
  zora:cutoff 1d-30
end

geometry noautosym noautoz
$(awk '$5="";FNR>2' $i )
end

task ${task}
EOF
done
done
