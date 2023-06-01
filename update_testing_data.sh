#!/bin/bash
ROOT=$(pwd)
TEST_CASES="test/test_cases.txt"

# TODO make dependent on set parameters and loop etc
[ -d "${ROOT}/inputs/gpref" ] && rm -rf ${ROOT}/inputs/gpref
[ -d "${ROOT}/inputs/abetaref" ] && rm -rf ${ROOT}/inputs/abetaref
[ -d "${ROOT}/inputs/ttrref" ] && rm -rf ${ROOT}/inputs/ttrref
cp -r test/test_inputs/* inputs

echo "Deleting current test_outputs"
rm -rf test/test_outputs/*

echo "Generating data for examples from ${ROOT}/${TEST_CASES}..."
# while read case; do echo $case ; done < ${ROOT}/${TEST_CASES}
while read case; do
    out_dir=$(python multiego.py $case | grep outputs/)
    mv $out_dir test/test_outputs
done < ${ROOT}/${TEST_CASES}
# python multiego.py < "$(cat ${ROOT}/${TEST_CASES} | xargs -0 -l -d \\n echo
# python multiego.py --system=gpref --egos rc > /dev/null
# python multiego.py --system=gpref --egos production --epsilon 0.35 --train_from md_ensemble > /dev/null
# python multiego.py --system=abetaref --egos production --epsilon 0.35 --train_from native_MD > /dev/null
# python multiego.py --system=ttrref --egos production --epsilon 0.225 --train_from native_MD fibril_MD --check_with fibril_check --inter_epsilon 0.3 > /dev/null

echo "Finished generating test outputs!"
