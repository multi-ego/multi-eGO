set -e
set -o pipefail

tar -zxf test_inputs/make_mat_ttr/hh.tgz -C test_inputs/make_mat_ttr/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_ttr/histo --target_top test_inputs/make_mat_ttr/topol_md.top --mego_top test_inputs/make_mat_ttr/topol_ref.top --cutoff 0.75 --inter --out  test_inputs/make_mat_ttr/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_ttr/histo --target_top test_inputs/make_mat_ttr/topol_md.top --mego_top test_inputs/make_mat_ttr/topol_ref.top --cutoff 0.75 --out  test_inputs/make_mat_ttr/
zdiff test_inputs/make_mat_ttr/intramat_1_1.ndx.gz test_outputs/make_mat_ttr/intramat_1_1.ndx.gz
zdiff test_inputs/make_mat_ttr/intermat_1_1.ndx.gz test_outputs/make_mat_ttr/intermat_1_1.ndx.gz

tar -zxf test_inputs/make_mat_popc/hh.tgz -C test_inputs/make_mat_popc/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat_popc/histo --target_top test_inputs/make_mat_popc/topol_md.top --mego_top test_inputs/make_mat_popc/topol_ref.top --cutoff 0.75 --out  test_inputs/make_mat_popc/
zdiff test_inputs/make_mat_popc/intramat_1_1.ndx.gz test_outputs/make_mat_popc/intramat_1_1.ndx.gz

