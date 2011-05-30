#!/bin/bash

# requires stdeb package from PyPI
python setup.py build sdist
py2dsc dist/Ctrax-`cat version.txt`.tar.gz
