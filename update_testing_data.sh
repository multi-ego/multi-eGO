#!/bin/bash

# TODO make dependent on set parameters and loop etc
cp -r test/test_inputs/* inputs

echo "Generating data for examples..."

python multiego.py --system=gpref --egos=rc > /dev/null
python multiego.py --system=gpref --egos=production --epsilon=0.35 --train_from=md_ensemble > /dev/null
python multiego.py --system=abetaref --egos=production --epsilon=0.35 --train_from=native_MD > /dev/null

rm -r test/test_outputs
cp -r outputs/gpref_rc outputs/gpref_production_e0.35_0.35 outputs/abeta_production_e0.35_0.35 test/test_outputs
