#!/bin/bash
set -eou pipefail

PARAMS=""
while (( "$#" )); do
  case "$1" in
    -d|--dev)
      shift 1
      DEVELOP=1
      ;;
    *) # preserve positional args
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
eval set -- "$PARAMS"

registry="registry.gitlab.com/shenglab/dose-calculation"
tag="stable"

fulltag="${registry}:${tag}"

echo "Building docker image..."
sudo docker build --tag "${fulltag}" .

echo "Pushing docker image to registry..."
sudo docker push "${fulltag}"

if ! [ ${DEVELOP:=} == 1 ]; then echo "Cleaning up development system"
  sudo docker image rm "${fulltag}"
fi
