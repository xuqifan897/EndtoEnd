#!/bin/bash

# launch debugger with args
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf
CUDA_DEBUGGER_SOFTWARE_PREEMPTION=1 cuda-gdb --args $PREPROCESS_FUNC_CALL $DEBUG_ARGS "$@"
