set -e
set -o pipefail
mkdir test_inputs/cmdata/histo
cd test_inputs/cmdata/histo
cmdata -f ../traj.xtc -s ../protein.tpr --mode=intra
cd ../../../
tar -zxf test_outputs/cmdata/hh.tgz -C test_outputs/cmdata/
for i in `ls -1 test_inputs/cmdata/histo`; do
  diff test_inputs/cmdata/histo/$i test_outputs/cmdata/histo/$i
done
