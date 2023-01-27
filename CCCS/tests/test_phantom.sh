#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf

PREPROCESS_FUNC_CALL="$PREPROCESS_BIN --density=$PHANTOM_VOL --fmaps=$PHANTOM_FMAPS --nphi=8 --ntheta=8 --voxsize=0.2 --ssfactor=1 --verbose --timing"
echo -e "$PREPROCESS_FUNC_CALL\n"
exec $PREPROCESS_FUNC_CALL
