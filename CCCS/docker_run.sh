#!/bin/bash

# Automatically run full-beam dose calculation and place result in /results
# To use, first build the docker image, then run it with no command, which will cause this script to be called
# The docker run command must contain two volume mappings for the input and output directories. Additionally,
#     the input directory must contain a './ct' folder and a './beamlist.txt' file which are already passed
#     to the command below

mkdir -p /data /results
cd /results
/dosecalc/bin/dosecalc-preprocess --dicom="/data/ct" --beamlist="/data/beamlist.txt" --config="/data/config.json" --verbose | tee /results/dosecalc-preprocess.log
/dosecalc/bin/dosecalc-beam --verbose | tee /results/dosecalc-beam.log
rm /results/Dose.h5
