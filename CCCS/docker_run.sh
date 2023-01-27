#!/bin/bash
set -eou pipefail

# Automatically run full-beam dose calculation and place result in /data/results
# To use, first build the docker image, then run it with no command, which will cause this script to be called
# The docker run command must contain two volume mappings for the input and output directories. Additionally,
#     the input directory must contain a './ct' folder and a './beamlist.txt' file which are already passed
#     to the command below

mkdir -p /data/results
cd /data/results
/dosecalc/bin/dosecalc-preprocess --dicom="/data/ct" --beamlist="/data/beamlist.txt" --config="/data/config.json" --verbose | tee dosecalc-preprocess.log
/dosecalc/bin/dosecalc-beamlet --verbose | tee dosecalc-beamlet.log
/dosecalc/bin/dosecalc-beam    --verbose | tee dosecalc-beam.log
rm -f /data/results/Dose.raw
