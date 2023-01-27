#!/bin/bash

# run full host-side memcheck using valgrind using full patient data
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf
valgrind --leak-check=yes $DOSE_BEAMLET_FUNC_CALL "$@" 2>&1 | tee valgrind_beamlet.log

