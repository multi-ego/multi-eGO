# check for black
pip show black
if [ $? -eq 0 ]; then
python -m black -t py39 -t py310 -t py311 -t py312 --line-length=127 .
else
echo "black module not found!"
echo "you can install it as pip install git+https://github.com/psf/black"
fi
flake8 . --count --extend-ignore E501,W605,E731,E402,E711 --show-source --statistics
if [ $? -eq 1 ]; then
flak="Fix the remaining formatting issue before pushing your code"
fi
echo "Running the regtests"
python test/run_tests.py >& /dev/null
if [ $? -eq 1 ]; then
  echo "Regtests failing, you should not push your code"
  echo $flak
else
  echo "Regtests passed"
  if [ -z "$flak" ]; then
    echo "Ready to push"
  else
    echo $flak
  fi
fi
