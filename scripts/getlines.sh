#!/bin/bash

# extracts the desired trinocular data for Chicago
# from synthetic data  files
#

cd /Users/rfabbri/lib/data/synthcurves-multiview-3d-dataset/spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object


frames="42 54 62"
#set -x
for i in $frames; do
  echo "// 2D points for frame $i"
  sed -n '621p;3012p;3390p' frame_00$i-pts-2d.txt
done

echo "// 3D points"
sed -n '621p;3012p;3390p' crv-3D-pts.txt

#for i in $frames; do
#  echo "// tangents for frame $i"
#  sed -n '3012p;3390p' frame_00$i-tgts-2d.txt
#done

for i in $frames; do
  echo "// extrinsics for frame $i"
  cat frame_00$i.extrinsic
done
