#!/bin/bash

# launch debugger with args
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf
CUDA_DEBUGGER_SOFTWARE_PREEMPTION=1 cuda-gdb --args $DOSE_BEAMLET_FUNC_CALL "$@"
