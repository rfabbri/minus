#!/bin/bash

# extracts the desired trinocular data for Chicago
# from synthetic data  files
#

cd /Users/rfabbri/lib/data/synthcurves-multiview-3d-dataset/spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object

echoerr() { echo "LOG $@" 1>&2; }
echoerr hello world

frames="42 54 62"
# off by one mistake by hongyi:
lines="3011 3389 620"
#set -x
for i in $frames; do
  echoerr "2D points for frame $i"
  for l in $lines; do
    sed -n "${l}p" frame_00$i-pts-2d.txt
  done
done

for i in $frames; do
  echoerr "2D tangents for frame $i"
  for l in $lines; do
    sed -n "${l}p" frame_00$i-tgts-2d.txt
  done
done

echoerr id0 id1
echo 0 1

echoerr "3D points"
for l in $lines; do
  sed -n "${l}p" crv-3D-pts.txt
done

#for i in $frames; do
#  echoerr "tangents for frame $i"
#  sed -n '3012p;3390p' frame_00$i-tgts-2d.txt
#done
echoerr "3D points"
head -n 2 calib.intrinsic

for i in $frames; do
  echoerr "extrinsics for frame $i"
  cat frame_00$i.extrinsic
done
