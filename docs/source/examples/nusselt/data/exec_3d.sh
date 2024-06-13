#!/bin/bash

set -x
set -e

lx=1.0e+0
ly=1.0e+0
lz=1.0e+0
Ra=5.0e+3
Pr=${1}

# set initial condition
cd initial_condition
make output
lx=${lx} \
ly=${ly} \
lz=${lz} \
glisize=16 \
gljsize=16 \
glksize=16 \
python main.py output
cd ..

# build and execute
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
Ra=${Ra} \
Pr=${Pr} \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts

ly=${ly} \
lz=${lz} \
Ra=${Ra} \
Pr=${Pr} \
python \
  docs/source/examples/nusselt/data/process.py \
  output/log/heat_transfer.dat \
  output/log/injected_squared_velocity.dat \
  output/log/dissipated_squared_velocity.dat \
  output/log/dissipated_squared_temperature.dat \
  artifacts/nusselt.png

echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci.txt
echo "Date :" $(date) >> artifacts/ci.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci.txt

