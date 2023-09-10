#!/bin/bash

set -x
set -e

Ra=${1}

# set initial condition
cd initial_condition
make output
lx=1.0e+0 \
ly=2.0e+0 \
glisize=128 \
gljsize=256 \
uniformx=false \
python main.py output
cd ..

# build and execute
sed -i "s/false/true/g" src/param/implicit.c && cat src/param/implicit.c
make all && make output
timemax=2.0e+2 \
wtimemax=3.0e+2 \
log_rate=5.0e-1 \
save_rate=1.0e+3 \
save_after=1.0e+3 \
stat_rate=1.0e+0 \
stat_after=1.0e+3 \
coef_dt_adv=0.95 \
coef_dt_dif=0.95 \
Ra=${Ra} \
Pr=7.0e+0 \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts
cp output/log/nusselt.dat artifacts/nusselt_${Ra}.dat
echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci_${Ra}.txt
echo "Date :" $(date) >> artifacts/ci_${Ra}.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci_${Ra}.txt

