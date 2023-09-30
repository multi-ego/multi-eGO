#!/bin/bash

i_block=$1
mat_1=$2
mat_2=$3

echo "Group $mat1 with $mat2 at atom $i_block"

paste $mat_1 $mat_2 | awk -v i=$i_block '{if(($2<=i && $4<=i) || ($2>i && $4>i)){print $1,$2,$3,$4,$5,$6,$7} else{print $8,$9,$10,$11,$12,$13,$14}}' > grouped_intramat.ndx

echo "Finished job"
