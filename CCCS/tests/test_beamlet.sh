#!/bin/bash

# run simple test using full patient data
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf
DOSE_BEAMLET_FUNC_CALL="$DOSE_BEAMLET_FUNC_CALL $@"

echo -e "$DOSE_BEAMLET_FUNC_CALL\n"
exec $DOSE_BEAMLET_FUNC_CALL
