#!/bin/bash
set -euo pipefail

# must match MIM structure names
BBox="BODY" # bounding box - dose is calculated in this box (speed up runtime)
dicomdata="dicomdata/"
configfile="config.json"
beamlist="beamlist.txt"
structures="structures.json"

# Quality Settings
voxelsize='0.25' # [units: cm]
sparsity='1e-4' # probably don't need to change ever

#-------------------------------------------------------------------------------------------------
export DOSECALC_DATA="/raid10/rs4pi/latest/data"

# call preprocess, save a log of the output automatically
dosecalc-preprocess \
    --dicom="$dicomdata" \
    --beamlist="$beamlist" \
    --structures="$structures" \
    --config="$configfile" \
    --bbox="$BBox" \
    --voxsize="$voxelsize" \
    2>&1 | tee "dosecalc-preprocess.log"

echo -e "\n\n=================================================================================\n\n"

# call dosecalc-beamlet, save a log of the output automatically
dosecalc-beamlet \
    --sparsify-threshold="$sparsity" \
    2>&1 | tee "dosecalc-beamlet.log"

echo -e "\n\n=================================================================================\n\n"	
