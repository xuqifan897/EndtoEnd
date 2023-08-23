#!/bin/bash

# resultFolder="/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14"

# if [ ! -d ${resultFolder} ]; then
#     mkdir ${resultFolder}
# fi

# logFile="${resultFolder}/log.txt"

# ( time ./build/tubScore \
#     --resultFolder /data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14 \
#     --nParticles 10000000 ) 2>&1 > ${logFile} &


if false; then
    # Now we are going to calculate the IPB kernels of different energy, 
    # to form a polychromatic kernel. In Ryan Neph's code, he implemented 
    # a 6MV photon beam. According to the spectrum file, the 6mv kernel 
    # contains the following monoenergetic beams:
    # 0.2 0.3 0.4 0.5 0.6 0.8 1.00 1.25 1.50 2.00 3.00 4.00 5.00 6.00

    motherFolder="/data/qifan/projects/EndtoEnd/results/spec6MV"
    if [ ! -d ${motherFolder} ];
    then
        mkdir ${motherFolder}
    fi

    for energy in 0.2 0.3 0.4 0.5 0.6 0.8 1.00 1.25 1.50 2.00 3.00 4.00 5.00 6.00;
    do
        resultFolder="${motherFolder}/E${energy}"
        if [ ! -d ${resultFolder} ];
        then
            mkdir ${resultFolder}
        fi

        logFile="${resultFolder}/log.txt"

        (time ./build/tubScore \
            --resultFolder ${resultFolder} \
            --Energy ${energy} \
            --nParticles 10000000 ) > ${logFile}
    done
fi

if false; then
    # This time, we are going to calculate water dose.
    motherFolder="/data/qifan/projects/EndtoEnd/results/waterIPB"
    if [ ! -d ${motherFolder} ]; then
        mkdir ${motherFolder}
    fi

    resultFolder="${motherFolder}/mono6MeV"
    if [ ! -d ${resultFolder} ]; then
        mkdir ${resultFolder}
    fi

    logFile="${resultFolder}/log.txt"

    (time ./build/tubScore \
        --resultFolder ${resultFolder} \
        --Energy 6.0 \
        --nParticles 10000000 ) > ${logFile}
fi

if true; then
    # This time, we are going to calculate bone dose
    motherFolder="/data/qifan/projects/EndtoEnd/results/waterIPB"
    if [ ! -d ${motherFolder} ]; then
        mkdir ${motherFolder}
    fi

    resultFolder="${motherFolder}/bone6MeV"
    if [ ! -d ${resultFolder} ]; then
        mkdir ${resultFolder}
    fi

    logFile="${resultFolder}/log.txt"

    (time ./build/tubScore \
        --resultFolder ${resultFolder} \
        --Energy 6.0 \
        --nParticles 10000000 ) > ${logFile}
fi