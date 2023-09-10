#!/bin/bash

set -x
set -e

Pr=${1}

# set initial condition
cd initial_condition
make output
lx=1.0e+0 \
ly=1.0e+0 \
lz=1.0e+0 \
glisize=16 \
gljsize=16 \
glksize=16 \
uniformx=false \
python main.py output
cd ..

# build and execute
sed -i "s/DNDIMS=2/DNDIMS=3/g" Makefile
sed -i "s/false/true/g" src/param/implicit.c
make all && make output
timemax=5.0e+3 \
wtimemax=6.0e+2 \
log_rate=1.0e+0 \
save_rate=5.0e+3 \
save_after=1.0e+3 \
stat_rate=1.0e+0 \
stat_after=1.0e+3 \
coef_dt_adv=0.95 \
coef_dt_dif=0.95 \
Ra=5.0e+3 \
Pr=${Pr} \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts
python \
  docs/source/examples/nu/data/process.py \
  output/log/nusselt.dat \
  artifacts/nusselt_${Pr}_3d.png
echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci_${Pr}_3d.txt
echo "Date :" $(date) >> artifacts/ci_${Pr}_3d.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci_${Pr}_3d.txt

