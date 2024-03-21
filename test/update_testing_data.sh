#!/bin/bash
MEGO_ROOT=..
TEST_DIR=${MEGO_ROOT}/test
TEST_CASES="test/test_cases.txt"

# TODO make dependent on set parameters and loop etc
[ -d "${MEGO_ROOT}/inputs/gpref" ] && rm -rf ${MEGO_ROOT}/inputs/gpref
[ -d "${MEGO_ROOT}/inputs/abetaref" ] && rm -rf ${MEGO_ROOT}/inputs/abetaref
[ -d "${MEGO_ROOT}/inputs/ttrref" ] && rm -rf ${MEGO_ROOT}/inputs/ttrref
[ -d "${MEGO_ROOT}/inputs/lyso-bnz_ref" ] && rm -rf ${MEGO_ROOT}/inputs/lyso-bnz_ref
cp -r test_inputs/* ${MEGO_ROOT}/inputs

echo "Deleting current test_outputs"
rm -rf test_outputs/gpref*
rm -rf test_outputs/abetaref*
rm -rf test_outputs/ttrref*
rm -rf test_outputs/lyso-bnz_ref*

echo "Generating data for examples from ${TEST_CASES}..."
cd .. 
# while read case; do echo $case ; done < ${ROOT}/${TEST_CASES}
while read case; do
    out_dir=$(python multiego.py $case | grep ${MEGO_ROOT}/outputs/)
    mv $out_dir test/test_outputs
done < ${TEST_CASES}

cd test 
echo "Finished generating test outputs!"
