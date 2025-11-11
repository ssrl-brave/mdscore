#!/bin/bash

prmtop=$1

cat > parmed.in << EOF
parm ${prmtop}
hmassrepartition
outparm system_hmass.top
EOF

parmed -i parmed.in | tee parmed.log
