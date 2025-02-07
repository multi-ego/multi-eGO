#!/bin/bash

set -e
set -o pipefail

python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_ttr/hh.tgz --target_top test_inputs/make_mat_ttr/topol_md.top --mego_top test_inputs/make_mat_ttr/topol_ref.top --cutoff 0.75 --mode intra+same --out  test_inputs/make_mat_ttr/ --tar
#diff <(gzip -dc test_inputs/make_mat_ttr/intramat_1_1.ndx.gz) <(gzip -dc test_outputs/make_mat_ttr/intramat_1_1.ndx.gz)
#diff <(gzip -dc test_inputs/make_mat_ttr/intermat_1_1.ndx.gz) <(gzip -dc test_outputs/make_mat_ttr/intermat_1_1.ndx.gz)

tar -zxf test_inputs/make_mat_popc/hh.tgz -C test_inputs/make_mat_popc/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_popc/histo --target_top test_inputs/make_mat_popc/topol_md.top --mego_top test_inputs/make_mat_popc/topol_ref.top --cutoff 0.75 --mode intra --out  test_inputs/make_mat_popc/ 
#diff <(gzip -dc test_inputs/make_mat_popc/intramat_1_1.ndx.gz) <(gzip -dc test_outputs/make_mat_popc/intramat_1_1.ndx.gz)

tar -zxf test_inputs/make_mat_het_trim/hh.tgz -C test_inputs/make_mat_het_trim/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_het_trim/histo --target_top test_inputs/make_mat_het_trim/topol_mego.top --mego_top test_inputs/make_mat_het_trim/topol.top --cutoff 0.75 --out  test_inputs/make_mat_het_trim/ 
#diff <(gzip -dc test_inputs/make_mat_het_trim/intramat_1_1.ndx.gz) <(gzip -dc test_outputs/make_mat_het_trim/intramat_1_1.ndx.gz)
#diff <(gzip -dc test_inputs/make_mat_het_trim/intramat_2_2.ndx.gz) <(gzip -dc test_outputs/make_mat_het_trim/intramat_2_2.ndx.gz)
#diff <(gzip -dc test_inputs/make_mat_het_trim/intramat_3_3.ndx.gz) <(gzip -dc test_outputs/make_mat_het_trim/intramat_3_3.ndx.gz)
#diff <(gzip -dc test_inputs/make_mat_het_trim/intermat_1_2.ndx.gz) <(gzip -dc test_outputs/make_mat_het_trim/intermat_1_2.ndx.gz)
#diff <(gzip -dc test_inputs/make_mat_het_trim/intermat_1_3.ndx.gz) <(gzip -dc test_outputs/make_mat_het_trim/intermat_1_3.ndx.gz)
#diff <(gzip -dc test_inputs/make_mat_het_trim/intermat_2_3.ndx.gz) <(gzip -dc test_outputs/make_mat_het_trim/intermat_2_3.ndx.gz)
