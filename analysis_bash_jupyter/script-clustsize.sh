#!/bin/bash

gmx_mpi clustsize -f prod_2000-1mM.xtc -s prod_2000-1mM.tpr -mol -hct -dt 1000000 -nice 0 -cut 0.4  &> clustlog.txt & 

exit