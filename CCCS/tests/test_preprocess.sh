#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/test_settings.conf
PREPROCESS_FUNC_CALL="$PREPROCESS_FUNC_CALL $@"

echo -e "$PREPROCESS_FUNC_CALL\n"
exec $PREPROCESS_FUNC_CALL
