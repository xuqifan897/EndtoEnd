#!/bin/bash

# run simple test using full patient data
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf
DOSE_BEAM_FUNC_CALL="$DOSE_BEAM_FUNC_CALL $@"

echo -e "$DOSE_BEAM_FUNC_CALL\n"
exec $DOSE_BEAM_FUNC_CALL
