#!/bin/bash

if (( $# != 1 )); then
    printf "usage: $0 <log-file>"  
    return 1
fi

echo "Avg Kernel Convolve Time:      " $( awk '/Kernel Convolve Time/ {IFS=" *"; print $7}' "$1" | awk 'BEGIN{sum=0} /\d*/ {sum+=$0} END{print sum/NR "ms"}' )
echo "Avg Write Sparse Beamlet Time: " $( awk '/Write Sparse Beamlet/ {IFS=" *"; print $9}' "$1" | awk 'BEGIN{sum=0} /\d*/ {sum+=$0} END{print sum/NR "s"}' )
