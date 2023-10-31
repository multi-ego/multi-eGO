#!/bin/bash
MEGO_ROOT=..
TEST_CASES="test_cases.txt"

# TODO make dependent on set parameters and loop etc
[ -d "${MEGO_ROOT}/inputs/gpref" ] && rm -rf ${MEGO_ROOT}/inputs/gpref
[ -d "${MEGO_ROOT}/inputs/abetaref" ] && rm -rf ${MEGO_ROOT}/inputs/abetaref
[ -d "${MEGO_ROOT}/inputs/ttrref" ] && rm -rf ${MEGO_ROOT}/inputs/ttrref
cp -r test_inputs/* ${MEGO_ROOT}/inputs

echo "Deleting current test_outputs"
rm -rf test_outputs/*

echo "Generating data for examples from ${TEST_CASES}..."
# while read case; do echo $case ; done < ${ROOT}/${TEST_CASES}
while read case; do
    out_dir=$(python ../multiego.py $case | grep ${MEGO_ROOT}/outputs/)
    mv $out_dir test_outputs
done < ${TEST_CASES}
# python multiego.py < "$(cat ${ROOT}/${TEST_CASES} | xargs -0 -l -d \\n echo
# python multiego.py --system=gpref --egos rc > /dev/null
# python multiego.py --system=gpref --egos production --epsilon 0.35 --train_from md_ensemble > /dev/null
# python multiego.py --system=abetaref --egos production --epsilon 0.35 --train_from native_MD > /dev/null
# python multiego.py --system=ttrref --egos production --epsilon 0.225 --train_from native_MD fibril_MD --check_with fibril_check --inter_epsilon 0.3 > /dev/null

echo "Finished generating test outputs!"
