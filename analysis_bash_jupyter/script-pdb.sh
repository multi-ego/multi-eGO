#!/bin/bash
# for-loopcmd.sh

# the two arrays are based on the maxclust.ndx -> check which clust is the one you want
array=("clust 1113", "clust 791", "clust 1199", "clust 1665", "clust 849", "clust 846", "clust 875", "clust 1035", "clust 1224", "clust 1390", "clust 1146", "clust 321", "clust 1194") # 33 is not fibril
array2=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "Z") # Z is not fibril


for ((i=0;i<${#array[@]}; ++i)); do
  echo -e ${array[i]}\n | gmx_mpi editconf -f prod_2000-1mM.tpr -n full-index.ndx -label ${array2[i]} -o clust${array[i]}.pdb
  echo -e ${array[i]}\n | gmx_mpi editconf -f prod_2000-1mM.tpr -n full-index.ndx -label ${array2[i]} -o clust${array[i]}.esp
  sed 's/{//' clust${array[i]}.esp | awk '{print $1+1}' > numero${array[i]}
  sed -i '1s/^/1\n/' numero${array[i]}
  paste numero${array[i]} clust${array[i]}.pdb > chain${array[i]}
  cat chain* | grep ATOM | sort -k1 -g > temp_label-chain.pdb
done

rm numero*
rm chain*
rm clust*
cut -f 1 --complement temp_label-chain.pdb > label-chain.pdb
mv label-chain.pdb chain-label.pdb

exit