##Workflow:

1) Use divide_intramat.sh (or domain_intramat.py) on the md intramat. This will generate an intramat (split_intramat.ndx) divided in blocks (0 for inter-domain contact) to be used for the inter-domain random coil
2) Use the rc intramat and the split_intramat.ndx to generate a mego ff with the intra-domain interaction on and the inter-domain interaction off. This will be the simulation of the inter-domain random coil.
3) Calculate the intramat of the inter-domain random coil
4) Use group_intramats.sh with rc intramat and inter-domain rc intramat. This will generate a new rc intramat with off diagonal blocks obtained from the inter-domain rc
5) Use the new rc_intramat and the full md intramat for multi-eGO to generate a new ff for multi-domain proteins


##Scrits:

divide_intramat.sh takes as input the atom defining the division between two domains and divide the second input (intramat) into four blocks setting to 0 all the inter-domain interaction

> bash divide_intramat.sh num(int) path_to_md_intramat

group_inrtamat.sh takes as input the atom defining the division between two domains, the path of the random_coil intramat and the path of the inter-domain random coil. The script regroups the intramat so that the intra-domain matrix blocks are the random_coil ones and the inter_domain matrix blocks are the inter-domain "random coil" ones

> bash group_intrmat.sh num(int) path_to_random_coil_intramat path_to_inter_domain_rc_intramat

Alternatively the python code (domain_intramat.py) can be used in the same way
