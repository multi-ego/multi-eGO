set -o pipefail
tar -zxvf test_inputs/make_mat/hh.tgz -C test_inputs/make_mat/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat/histo --target_top test_inputs/make_mat/topol_md.top --mego_top test_inputs/make_mat/topol_ref.top --cutoff 0.75 --inter --out  test_inputs/make_mat/
python ../tools/make_mat/make_mat.py --histo test_inputs/make_mat/histo --target_top test_inputs/make_mat/topol_md.top --mego_top test_inputs/make_mat/topol_ref.top --cutoff 0.75 --out  test_inputs/make_mat/
diff test_inputs/make_mat/intramat_1_1.ndx test_outputs/make_mat/intramat_1_1.ndx
diff test_inputs/make_mat/intermat_1_1.ndx test_outputs/make_mat/intermat_1_1.ndx

