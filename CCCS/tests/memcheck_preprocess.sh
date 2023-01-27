#!/bin/bash

# run full memcheck using full patient data
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf
cuda-memcheck --leak-check full --report-api-errors all --racecheck-report all $PREPROCESS_FUNC_CALL "$@" | tee memcheck_preprocess.log
