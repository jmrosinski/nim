#!/bin/bash
#JR -x means debug mode
set -x
tail=`basename "$PWD"`
if [[ $tail != "run" || ! -d $PWD ]]; then
   echo "You must be in the run/ directory to run $0"
   exit 1
fi

nimrelpath=../src/dynamics/nim
if [ test ! -x $nimrelpath ]; then
  echo "nim either not found or not executable"
  exit 1
fi

rundir=../run_$$
if [ -d $rundir ]; then
  echo "directory $rundir already exists"
  exit 1
fi

mkdir $rundir
cp -p $nimrelpath ./NIMnamelist $rundir
cd $rundir
ln -s /data1/spiralOverlap/data/FV3DAT/G3D_G05K096.dat ./g3d.dat
ln -s /data1/spiralOverlap/data/FV3DAT/STD_G05K096.txt ./STD.txt
ln -s /data1/spiralOverlap/data/ATMDAT/ini_G05K096.dat ./ini.dat
./nim
exit 0
