#!/bin/bash

i_block=$1
mat=$2

echo "Dividing $2 at atom i = $i_block"
echo 

awk -v i=$i_block '{if(($2<=i && $4<=i) || ($2>i && $4>i)){print $0} else{print $1,$2,$3,$4,"0.000000","0.000000+e00" ,$7}}' $2 > split_intramat.ndx

echo "Finished job"
