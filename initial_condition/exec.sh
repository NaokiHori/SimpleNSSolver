#!/bin/bash

## domain
# domain lengths
export lx=1.0e+0
export ly=2.0e+0
# number of cell centers
export glisize=128
export gljsize=256
# grid is uniformly distributed in x (true) or not (false)
export uniformx=false

## where to write resulting NPY files
dirname="output"

python3 main.py ${dirname}
