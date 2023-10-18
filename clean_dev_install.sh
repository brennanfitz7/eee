#!/bin/bash

find . | grep -E "(/__pycache__$|\.pyc$|\.pyo$)" | xargs rm -rf
python setup.py develop
