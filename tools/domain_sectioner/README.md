divide_inrtamat.sh takes as input the atom defining the division between two domains and divide the second input (intramat) into four blocks setting to 0 all the inter-domain interaction

> bash divide_intramat.sh num(int) path_to_md_intramat

group_inrtamat.sh takes as input the atom defining the division between two domains, the path of the random_coil intramat and the path of the inter-domain random coil. The script regroups the intramat so that the intra-domain matrix blocks are the random_coil ones and the inter_domain matrix blocks are the inter-domain "random coil" ones

> bash group_intrmat.sh num(int) path_to_random_coil_intramat path_to_inter_domain_rc_intramat

Alternatively the python code can be used in the same way
