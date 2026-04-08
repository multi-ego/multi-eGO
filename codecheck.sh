# check for black
pip show black
if [ $? -eq 0 ]; then
  python -m black -t py310 -t py311 -t py312 -t py313 --line-length=120 .
else
  echo "black module not found!"
  echo "you can install it as pip install git+https://github.com/psf/black"
fi
flake8 . --count --extend-ignore E501,W605,E731,E402,E711 --show-source --statistics
if [ $? -eq 1 ]; then
  flak="Fix the remaining formatting issue before pushing your code"
  echo $flak
  exit 1
fi
echo "Running unit tests part 1"
pytest tests/test_multiego.py
if [ $? -eq 1 ]; then
  echo "UnitTests failing, you should not push your code"
  exit 1
else
  echo "UnitTest passed"
fi
echo "Running unit tests part 2"
pytest tests/test_apply_symmetries.py
pytest tests/test_symmetrize_reference_contacts.py
if [ $? -eq 1 ]; then
  echo "UnitTests failing, you should not push your code"
  exit 1
else
  echo "UnitTest passed"
fi
echo "Running the regtests: mego"
pytest tests/run_tests.py
if [ $? -eq 1 ]; then
  echo "Regtests failing, you should not push your code"
  exit 1
else
  echo "Regtests passed"
fi
echo "Running the regtests: make_mat"
pytest tests/test_make_mat.py
if [ $? -eq 1 ]; then
  echo "Regtests failing, you should not push your code"
  exit 1
else
  echo "Regtests passed"
fi
echo "All tests passed. Ready to push"
