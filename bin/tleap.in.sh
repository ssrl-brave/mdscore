#!/bin/bash

LIG=$1
FFMOD=$2
RECEPTOR=$3

cat > tleap.in << EOF
source oldff/leaprc.ff99SB
source leaprc.gaff2
UNL=loadmol2 ${LIG}
loadamberparams ${FFMOD}
saveoff UNL ligand.lib
saveamberparm UNL ligand.prmtop ligand.rst7
quit
EOF
tleap -f tleap.in | tee tleap.log

cat > tleap_0.in << EOF
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
prot = loadpdb ${RECEPTOR}
loadamberparams ${FFMOD}
lig=loadmol2 ${LIG}
loadoff ligand.lib
complex=combine {prot lig}
addions complex Na+ 0
addions complex Cl- 0
solvatebox complex OPCBOX 8.0
savepdb complex complex_solvated.pdb
saveamberparm complex complex_solvated.prmtop complex_solvated.inpcrd
quit
EOF

tleap -f tleap_0.in | tee tleap_0.log
