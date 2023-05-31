#!/bin/bash
ROOT=$(pwd)

# TODO make dependent on set parameters and loop etc
[ -d "${ROOT}/inputs/gpref" ] && rm -rf ${ROOT}/inputs/gpref
[ -d "${ROOT}/inputs/abetaref" ] && rm -rf ${ROOT}/inputs/abetaref
[ -d "${ROOT}/inputs/ttrref" ] && rm -rf ${ROOT}/inputs/ttrref
cp -r test/test_inputs/* inputs

echo "Generating data for examples..."
python multiego.py --system=gpref --egos rc > /dev/null
python multiego.py --system=gpref --egos production --epsilon 0.35 --train_from md_ensemble > /dev/null
python multiego.py --system=abetaref --egos production --epsilon 0.35 --train_from native_MD > /dev/null
python multiego.py --system=ttrref --egos production --epsilon 0.225 --train_from native_MD fibril_MD --check_with fibril_check --inter_epsilon 0.3 > /dev/null

echo "Finished generating\nMoving data to reference output directory"

rm -rf test/test_outputs/*
mv outputs/gpref_rc outputs/gpref_production_e0.35_0.35 outputs/abetaref_production_e0.35_0.35 outputs/ttrref_production_e0.225_0.3 test/test_outputs
