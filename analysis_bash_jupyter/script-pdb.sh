#!/bin/bash
# for-loopcmd.sh 

# the two arrays are based on the maxclust.ndx -> check which clust is the one you want
array=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "43")
array2=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "Z") 
 

for ((i=0;i<${#array[@]}; ++i)); do
	 echo -e ${array[i]} \n | gmx_mpi editconf -f prod_2000-1mM.tpr -n pdb.ndx -label ${array2[i]} -o clust${array[i]}.pdb
	 echo -e ${array[i]} \n | gmx_mpi editconf -f prod_2000-1mM.tpr -n pdb.ndx -label ${array2[i]} -o clust${array[i]}.esp
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

temp_label-chain.pdb
exit
