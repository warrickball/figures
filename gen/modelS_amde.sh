#!/usr/bin/env bash

# bash script to generate Model S eigenfunction data
# assumes that ADIPLS executables are in $PATH

if [ ! -f ../data/modelS.fgong ]
then
    curl http://users-phys.au.dk/jcd/solar_models/fgong.l5bi.d.15c -o ../data/modelS.fgong
fi
fgong-amdl.d ../data/modelS.fgong modelS.amdl
redistrb.c.d redistrb.c.in
adipls.c.d adipls.c.in
mv -v modelS.r.amdl modelS.amde ../data/
rm adipls-status.log fort.9 ttt.* modelS.amdl 
