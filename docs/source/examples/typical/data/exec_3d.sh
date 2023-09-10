#!/bin/bash

set -x
set -e

# set initial condition
cd initial_condition
make output
lx=1.0e+0 \
ly=1.0e+0 \
lz=1.0e+0 \
glisize=32 \
gljsize=32 \
glksize=32 \
uniformx=false \
python main.py output
cd ..

# build and execute
sed -i "s/DNDIMS=2/DNDIMS=3/g" Makefile
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
Ra=1.0e+6 \
Pr=1.0e+1 \
mpirun -n 2 --oversubscribe ./a.out initial_condition/output

# post process
mkdir artifacts
python \
  docs/source/examples/typical/data/snapshot_3d.py \
  $(find output/save -type d | sort | tail -n 1) \
  artifacts/snapshot_3d.png
python \
  docs/source/examples/typical/data/divergence.py \
  output/log/divergence.dat \
  artifacts/divergence_3d.png
python \
  docs/source/examples/typical/data/nusselt_time.py \
  output/log/nusselt.dat \
  artifacts/nusselt_time_3d.png
python \
  docs/source/examples/typical/data/nusselt_x.py \
  $(find output/stat -type d | sort | tail -n 1) \
  artifacts/nusselt_x_3d.png
python \
  docs/source/examples/typical/data/std.py \
  $(find output/stat -type d | sort | tail -n 1) \
  artifacts/std_3d.png
echo "OS   :" $(cat /etc/os-release | grep PRETTY_NAME | awk -F "=" '{print $2}') >> artifacts/ci_3d.txt
echo "Date :" $(date) >> artifacts/ci_3d.txt
echo "Hash :" $(git rev-parse HEAD) >> artifacts/ci_3d.txt

