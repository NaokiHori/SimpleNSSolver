#!/bin/bash

set -x
set -e

lx=1.0e+0
ly=2.0e+0
Ra=1.0e+8
Pr=1.0e+1

# set initial condition
cd initial_condition
make output
lx=${lx} \
ly=${ly} \
glisize=128 \
gljsize=256 \
python main.py output
cd ..

# build and execute
make all && make output
timemax=3.0e+2 \
wtimemax=6.0e+2 \
log_rate=5.0e-1 \
save_rate=5.0e+1 \
save_after=0.0e+0 \
stat_rate=1.0e-1 \
stat_after=2.0e+2 \
coef_dt_adv=0.95 \
coef_dt_dif=0.95 \
Ra=${Ra} \
Pr=${Pr} \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts

# visualise
python \
  docs/source/examples/typical/data/snapshot_2d.py \
  $(find output/save -type d | sort | tail -n 1) \
  artifacts/snapshot.png

# divergence-t plot
python \
  docs/source/examples/typical/data/divergence.py \
  output/log/max_divergence.dat \
  artifacts/divergence.png

# Nu-t plot
ly=${ly} \
Ra=${Ra} \
Pr=${Pr} \
python \
  docs/source/examples/typical/data/nusselt_time.py \
  output/log/heat_transfer.dat \
  output/log/injected_squared_velocity.dat \
  output/log/dissipated_squared_velocity.dat \
  output/log/dissipated_squared_temperature.dat \
  artifacts/nusselt_time.png

# Nu-x plot
ly=${ly} \
Ra=${Ra} \
Pr=${Pr} \
python \
  docs/source/examples/typical/data/nusselt_x.py \
  $(find output/stat -type d | sort | tail -n 1) \
  artifacts/nusselt_x.png

# standard deviations in x
python \
  docs/source/examples/typical/data/std.py \
  $(find output/stat -type d | sort | tail -n 1) \
  artifacts/std.png

echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci.txt
echo "Date :" $(date) >> artifacts/ci.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci.txt

