#!/bin/bash

set -x
set -e

lx=1.0e+0
ly=2.0e+0
Ra=1.0e+4
Pr=${1}

# set initial condition
cd initial_condition
make output
lx=${lx} \
ly=${ly} \
glisize=32 \
gljsize=64 \
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
Ra=${Ra} \
Pr=${Pr} \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts

ly=${ly} \
Ra=${Ra} \
Pr=${Pr} \
python \
  docs/source/examples/nusselt/data/process.py \
  output/log/heat_transfer.dat \
  output/log/injected_squared_velocity.dat \
  output/log/dissipated_squared_velocity.dat \
  output/log/dissipated_squared_temperature.dat \
  artifacts/nusselt.png

echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci_${Pr}.txt
echo "Date :" $(date) >> artifacts/ci_${Pr}.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci_${Pr}.txt

