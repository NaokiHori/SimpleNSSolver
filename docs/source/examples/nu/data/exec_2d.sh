#!/bin/bash

set -x
set -e

Pr=${1}

# set initial condition
cd initial_condition
make output
lx=1.0e+0 \
ly=2.0e+0 \
glisize=32 \
gljsize=64 \
uniformx=false \
python main.py output
cd ..

# build and execute
make all && make output
timemax=1.0e+3 \
wtimemax=6.0e+2 \
log_rate=1.0e+0 \
save_rate=1.0e+3 \
save_after=1.0e+3 \
stat_rate=1.0e+0 \
stat_after=1.0e+3 \
coef_dt_adv=0.95 \
coef_dt_dif=0.95 \
Ra=1.0e+4 \
Pr=${Pr} \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts
python \
  docs/source/examples/nu/data/process.py \
  output/log/nusselt.dat \
  artifacts/nusselt_${Pr}_2d.png
echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci_${Pr}_2d.txt
echo "Date :" $(date) >> artifacts/ci_${Pr}_2d.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci_${Pr}_2d.txt

