#!/bin/bash
MEGO_ROOT=$(realpath ..)
TEST_DIR=${MEGO_ROOT}/test
TEST_CASES="${TEST_DIR}/test_cases.txt"

# TODO make dependent on set parameters and loop etc
echo "Deleting current inputs"
[ -d "${MEGO_ROOT}/inputs/gpref" ] && rm -rf ${MEGO_ROOT}/inputs/gpref
[ -d "${MEGO_ROOT}/inputs/abetaref" ] && rm -rf ${MEGO_ROOT}/inputs/abetaref
[ -d "${MEGO_ROOT}/inputs/ttrref" ] && rm -rf ${MEGO_ROOT}/inputs/ttrref
[ -d "${MEGO_ROOT}/inputs/lyso-bnz_ref" ] && rm -rf ${MEGO_ROOT}/inputs/lyso-bnz_ref
cp -r test_inputs/* ${MEGO_ROOT}/inputs

# delete output directories
echo "Deleting current outputs"
[ -d "${MEGO_ROOT}/outputs/gpref" ] && rm -rf ${MEGO_ROOT}/outputs/gpref
[ -d "${MEGO_ROOT}/outputs/abetaref" ] && rm -rf ${MEGO_ROOT}/outputs/abetaref
[ -d "${MEGO_ROOT}/outputs/ttrref" ] && rm -rf ${MEGO_ROOT}/outputs/ttrref
[ -d "${MEGO_ROOT}/outputs/lyso-bnz_ref" ] && rm -rf ${MEGO_ROOT}/outputs/lyso-bnz_ref

echo "Deleting current test_outputs"
rm -rf ${TEST_DIR}/test_outputs/gpref
rm -rf ${TEST_DIR}/test_outputs/abetaref
rm -rf ${TEST_DIR}/test_outputs/ttrref
rm -rf ${TEST_DIR}/test_outputs/lyso-bnz_ref

echo "Generating data for examples from ${TEST_CASES}..."
while read case; do
    out_dir=$(python ${MEGO_ROOT}/multiego.py $case | grep "Output files written in" | awk '{print $5}')
    out_dir=$(realpath $out_dir)
    system=$(dirname $out_dir)
    system=$(basename $system)
    echo "Copying $out_dir to ${TEST_DIR}/test_outputs/${system}"
    mkdir -p ${TEST_DIR}/test_outputs/${system}
    cp -r $out_dir ${TEST_DIR}/test_outputs/${system}
done < ${TEST_CASES}

echo "Finished generating test outputs!"
