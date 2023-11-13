#!/bin/bash

for struc in $(ls *pdb)
do
    sed -i "s/\(.................POPC\)./\1A/g" ${struc}
    sed -i "s/\(.................POPS\)./\1B/g" ${struc}
    sed -i "s/\(.................POP2\)./\1C/g" ${struc}
    sed -i "s/\(.................CHOL\)./\1D/g" ${struc}
    sed -i "s/\(................. [A-Z][A-Z][A-Z]\)./\1E/g" ${struc}
done
