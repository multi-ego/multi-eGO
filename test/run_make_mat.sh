#!/bin/bash

set -e
set -o pipefail

tar -zxf test_inputs/make_mat_ttr/hh.tgz -C test_inputs/make_mat_ttr/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_ttr/histo --target_top test_inputs/make_mat_ttr/topol_mego.top --mego_top test_inputs/make_mat_ttr/topol.top --cutoff 0.75 --mode intra+same --out  test_inputs/make_mat_ttr/
python ../tools/make_mat/HDF52ndx.py --input test_inputs/make_mat_ttr/intramat_1_1.ndx.h5 
python ../tools/make_mat/HDF52ndx.py --input test_inputs/make_mat_ttr/intermat_1_1.ndx.h5 
diff test_inputs/make_mat_ttr/intramat_1_1.ndx test_outputs/make_mat_ttr/intramat_1_1.ndx
diff test_inputs/make_mat_ttr/intermat_1_1.ndx test_outputs/make_mat_ttr/intermat_1_1.ndx

tar -zxf test_inputs/make_mat_popc/hh.tgz -C test_inputs/make_mat_popc/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_popc/histo --target_top test_inputs/make_mat_popc/topol_md.top --mego_top test_inputs/make_mat_popc/topol_ref.top --cutoff 0.75 --mode intra --out  test_inputs/make_mat_popc/ --noh5 
#diff <(gzip -dc test_inputs/make_mat_popc/intramat_1_1.ndx.h5) <(gzip -dc test_outputs/make_mat_popc/intramat_1_1.ndx.h5)
