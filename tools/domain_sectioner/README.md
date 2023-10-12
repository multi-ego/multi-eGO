##Workflow:

1) Use multi_domain_intramats.py with "--type split" on the md intramat. This will generate an intramat (split_intramat.ndx) divided in blocks (0 for inter-domain contact) to be used for the inter-domain random coil
2) Use the rc intramat and the split_intramat.ndx to generate a mego ff with the intra-domain interaction on and the inter-domain interaction off. This will be the simulation of the inter-domain random coil.
3) Calculate the intramat of the inter-domain random coil
4) Use multi_domain_intramats.py with "--type group" on rc intramat and inter-domain rc intramat. This will generate a new rc intramat with off diagonal blocks obtained from the inter-domain rc
5) Use the new rc_intramat and the full md intramat for multi-eGO to generate a new ff for multi-domain proteins


##Scrits:

type split:

> multi_domain_intramats.py --type split --md_intra PATH/md_intramat --target_top PATH/target_top --mego_top PATH/mego_top --out ouput_directory --dom_res idx,idx,...

This will split md intramat at the atoms corresponding to the beginning of the residues passed with --dom_res.


type group:

> multi_domain_intramats.py --type group --rc_intra PATH/rc_intramat --dom_rc_intra PATH/dom_rc_intramat --target_top PATH/target_top --mego_top PATH/mego_top --out ouput_directory --dom_res idx,idx,...

This will group rc intramat and domain rc intramat at the atoms corresponding to the beginning of the residues passed with --dom_res.