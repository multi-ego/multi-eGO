set -e
set -o pipefail

# Relative and absolute tolerances for numerical comparison.
# Adjust RTOL / ATOL here to tighten or loosen the test.
RTOL=${CMDATA_RTOL:-5e-4}
ATOL=${CMDATA_ATOL:-1e-7}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
COMPARE="$SCRIPT_DIR/compare_histo.py"

mkdir -p test_inputs/cmdata/histo
cd test_inputs/cmdata/histo
OMP_NUM_THREADS=2 cmdata -f ../traj.xtc -s ../protein.tpr --noh5
cd ../../../

tar -zxf test_outputs/cmdata/hh.tgz -C test_outputs/cmdata/

failures=0
for i in $(ls -1 test_inputs/cmdata/histo); do
  if ! python3 "$COMPARE" \
        test_inputs/cmdata/histo/$i \
        test_outputs/cmdata/histo/$i \
        --rtol "$RTOL" --atol "$ATOL"; then
    failures=$((failures + 1))
  fi
done

if [ "$failures" -gt 0 ]; then
  echo "FAILED: $failures file(s) exceeded tolerance (rtol=$RTOL, atol=$ATOL)"
  exit 1
fi
echo "PASSED: all files agree within tolerance (rtol=$RTOL, atol=$ATOL)"
