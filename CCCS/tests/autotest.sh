#!/bin/bash
#
# FOR FULL LIST OF DEPENCENCIES, SEE README.md
#
# Test script for validating data directory structure, build
# procedure, and execution of GPU based preprocessing and dose
# calculation code.
#
set -e -o pipefail

test_build_cmake() {
    local b=build
    mkdir -p "$b"
    test_run "$b" cmake ..
    test_run "$b" make
}

test_existence_file() {
    if [[ ! -f $1 ]]; then
        echo "$1: does not exist"
        return 1
    fi
}

test_run() {
    local dir=$1
    shift
    local op=( "$@" )
    local cmd=${op[0]}
    if [[ $cmd =~ \.sh$ ]]; then
        op=( bash "${op[@]}" )
    fi
    local log="$dir/$cmd.log"
    echo -n "Testing $cmd....."
    if ! (cd "$dir" && "${op[@]}") >& "$log"; then
        echo failure
        cat "$log"
        return 1
    fi
    echo success
}

test_main() {
    cd "$(dirname "${BASH_SOURCE[0]}")"/..
    mkdir -p data/{ctdata,dsa,images,results,spectra,temp}
    #setup dosecalc_testdata001_prostate setup needed

    echo FILE EXISTENCE TESTS:
    local test_files=(
        data/dsa/scaf6000
        data/dsa/radii.dat
        data/dsa/polar.dat
        data/spectra/spec_6mv.spec
        'data/ctdata/dosecalc_testdata001_prostate/2.16.840.1.114362.1.6.6.5.16628.11442066728.464245670.819.211.dcm'
        tests/input_files/beamlist_4pi.txt
        tests/input_files/structures_prostate.json
    )
    for f in ${test_files[@]}; do
        test_existence_file "$f"
    done
    echo success
    echo ''

    echo BUILD TESTS:
    test_build_cmake
    echo ''

    echo EXECUTION TESTS:
    test_run . tests/test_preprocess.sh
    test_run . tests/test_beamlet.sh
    test_run . tests/test_beam.sh
    echo All Tests Passed
}

test_main
