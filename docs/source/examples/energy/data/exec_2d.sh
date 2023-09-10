#!/bin/bash

set -x
set -e

# set initial condition
cd initial_condition
make output
lx=1.0e+0 \
ly=1.0e+0 \
glisize=32 \
gljsize=32 \
uniformx=false \
python main.py output
cd ..

# configure
export wtimemax=3.0e+2
export save_rate=1.0e+3
export save_after=1.0e+3
export stat_rate=1.0e+0
export stat_after=1.0e+3
export coef_dt_dif=0.95
export Ra=1.0e+100
export Pr=1.0e+0

# create random fields
make all && make output
timemax=5.0e+1 \
log_rate=5.0e-1 \
coef_dt_adv=0.95 \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# stash last flow field
dirname_ic=$(find output/save -type d | sort | tail -n 1)
mv ${dirname_ic} dirname_ic

# check decays for different time steps
for factor_adv in 0.1 0.2 0.4 0.8
do
  sed -i "s/param_add_buoyancy.*$/param_add_buoyancy = false;/g" src/param/buoyancy.c && cat src/param/buoyancy.c
  make datadel && make all && make output
  timemax=6.0e+1 \
  log_rate=1.0e-1 \
  coef_dt_adv=${factor_adv} \
  mpirun -n 2 --oversubscribe ./a.out dirname_ic
  mv output/log/energy.dat ./energy-${factor_adv}.dat
done

# post process
mkdir artifacts
python \
  docs/source/examples/energy/data/process.py \
  . \
  artifacts/energy1_2d.png \
  artifacts/energy2_2d.png
echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci_2d.txt
echo "Date :" $(date) >> artifacts/ci_2d.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci_2d.txt

