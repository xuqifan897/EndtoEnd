"""
This python script compares the number of particle 
entries and number of tracks. We assume they are equal
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def PartiTrajComp():
    RE01path = '/home/qifan/projects/EndtoEnd4/geant4/RE01'
    myOutput = os.path.join(RE01path, 'myOutput.txt')
    with open(myOutput, 'r') as f:
        lines = f.readlines()

    partiLineStart = 1382
    partiLineEnd = 1793
    trajLineStart = 1796
    trajLineEnd = 10264

    partiLines = lines[partiLineStart:partiLineEnd]
    trajLines = lines[trajLineStart:trajLineEnd]

    # # have a check
    # print(partiLines[0])
    # print(partiLines[-1])
    # print(trajLines[0])
    # print(trajLines[-1])

    # build particle records (trackID, parentID)
    partiRecords = []
    hierarchyStack = [0]
    for line in partiLines:
        # firstly, count the number of leading spaces
        _line = line.split('>')[1]
        _line = _line.split('==')[0]
        nSpaces = len(_line)
        hierarchy = int(nSpaces / 3)

        # secondly, extract trackID
        trackID = int(line.split(' ')[-1])
        hierarchyStack = hierarchyStack[:hierarchy] + [trackID]
        partiRecords.append((hierarchyStack[-1], hierarchyStack[-2]))
    # print(partiRecords)

    # build track records (trackID, parentID)
    trajRecords = []
    for line in trajLines:
        if 'TrackID' in line:
            parts = line.split('=')
            TrackIDPart = parts[1]
            ParentIDPart = parts[2]
            TrackID = int(TrackIDPart.split(':')[0])
            ParentID = int(ParentIDPart.split(':')[0])
            trajRecords.append((TrackID, ParentID))
    # print(trajRecords)

    # compare the trackIDs
    partiTracks = [a[0] for a in partiRecords]
    trajTracks = [a[0] for a in trajRecords]

    # firstly, see if every partiTrack is in trajTrack
    # the answer is false, in fact, there is an overlap
    # between partiTracks and trajTracks
    falseList = []
    for a in partiTracks:
        if a not in trajTracks:
            falseList.append(a)
    # print(falseList)
    # print(len(partiTracks), len(trajTracks), len(falseList))

    # secondly, check if trajRecords is a tree
    # the answer is false, in fact, there are some nodes whose 
    # parent node does not exist
    trajRecordsCopy = trajRecords.copy()
    newNodes = [0]
    while(True):
        newNodesNew = []
        for newNode in newNodes:
            # filter out the entries whose parentID == newNode
            filtered = [a for a in trajRecordsCopy if a[1] == newNode]
            # add the filtered trackID to newNodesNew
            newNodesNew.extend([a[0] for a in filtered])
            trajRecordsCopy = [a for a in trajRecordsCopy if a not in filtered]
        newNodes = newNodesNew
        if len(newNodes) == 0:
            break
    # select the parentIDs of the remaining trajRecordsCopy
    remaining = [a[1] for a in trajRecordsCopy]
    print(remaining)


def waterKernelOutputProcessing():
    """
    This function processes the output by the PHASE 0 of the code.
    The code runs 32 particles, which are distributed to several threads.
    """
    outputFolder = '/data/qifan/projects/EndtoEnd4/results/point6MeV32'
    outputFile = os.path.join(outputFolder, 'myOutput.txt')
    with open(outputFile, 'r') as f:
        lines = f.readlines()
    logStartLine = 974
    logEndLine = 3371
    logLines = lines[logStartLine:logEndLine]
    groups = {}
    for line in logLines:
        parts = line.split(' ')
        threadID = parts[0]
        if threadID not in groups:
            groups[threadID] = []
        groups[threadID].append(line)
    for key, item in groups.items():
        threadFile = os.path.join(outputFolder, key+'.txt')
        with open(threadFile, 'w') as f:
            f.writelines(item)


def examineIPBdose():
    """
    This function examines the output file of the IPB dose.
    """
    resultsFolder = '/data/qifan/projects/EndtoEnd4/results/figures'
    if not os.path.isdir(resultsFolder):
        os.mkdir(resultsFolder)

    flag = 2
    if flag == 0:
        # process point dose kernel
        binaryFile = '/data/qifan/projects/EndtoEnd4/results' \
            '/point6MeV1e7Rand/array.bin'
        shape = (199, 199, 400)
        array = np.fromfile(binaryFile, dtype=np.double)
        array = np.reshape(array, shape)

        # get the depth dose profile
        depthDose = np.sum(array, axis=(0, 1))
        xAxis = np.arange(-50, 350) * 0.1  # cm
        plt.plot(xAxis, depthDose)
        plt.xlabel('depth w.r.t. the first interaction point [cm]')
        plt.ylabel('energy deposition (a.u.)')
        plt.title('point kernel plot')
        outputFile = os.path.join(resultsFolder, 'pointKernel.png')
        plt.savefig(outputFile)
        plt.clf()
    elif flag == 1:
        # process IPB kernel
        binaryFile = '/data/qifan/projects/EndtoEnd4/results' \
            '/IPB6MeV1e7/array.bin'
        shape = (199, 199, 400)
        array = np.fromfile(binaryFile, dtype=np.double)
        array = np.reshape(array, shape)

        depthDose = np.sum(array, axis=(0, 1))
        xAxis = np.arange(400) * 0.1  # cm
        plt.plot(xAxis, depthDose)
        plt.xlabel('depth w.r.t. water surface [cm]')
        plt.ylabel('energy deposition (a.u.)')
        plt.title('infinitesimal pencil beam kernel plot')
        outputFile = os.path.join(resultsFolder, 'IPBKernel.png')
        plt.savefig(outputFile)
        plt.clf()
    elif flag == 2:
        pointKernelBinary = '/data/qifan/projects/EndtoEnd4/results' \
            '/point6MeV1e7Rand/array.bin'
        shape = (199, 199, 400)
        array = np.fromfile(pointKernelBinary, dtype=np.double)
        array = np.reshape(array, shape)
        reducedPointKernel = np.sum(array, axis=(0, 1))

        # according to Attix, the mass attenuation coefficient of
        # 6 MeV photon in water is 0.0277 m^2/kg. We assume that 
        # the density of water is 1 kg/m^3. The mass energy transfer
        # coefficient of 6 MeV photon in water is 0.0185 m^2/kg

        muOverRho = 0.0277  # cm^2/g
        # muOverRho = 0.0185  # cm^2/g
        mu = muOverRho * 0.1 * 1000
        depth = np.arange(400) * 1e-3  # m
        Terma = np.exp(-depth * mu)
        # plt.plot(depth, Terma)
        # plt.show()

        # convolve point kernel with Terma
        kernelSize = 400
        kernelBase = 50
        kernelMargin = kernelSize - kernelBase
        doseSize = 400
        kernelMargin = kernelSize - kernelBase
        dose = np.zeros(kernelSize)
        for i in range(doseSize):
            kernelIdxBegin = np.max((0, kernelBase - i))
            kernelIdxEnd = np.min((kernelSize, doseSize - i + kernelBase))
            
            doseIdxBegin = np.max((0, i-kernelBase))
            doseIdxEnd = np.min((doseSize, i + kernelMargin))

            assert kernelIdxEnd - kernelIdxBegin == doseIdxEnd - doseIdxBegin
            dose[doseIdxBegin : doseIdxEnd] += Terma[i] * \
                reducedPointKernel[kernelIdxBegin : kernelIdxEnd]
        
        # compare with the experimental results
        IPBKernelBinary = '/data/qifan/projects/EndtoEnd4/results/' \
            'IPB6MeV1e7/array.bin'
        array = np.fromfile(IPBKernelBinary, dtype=np.double)
        array = np.reshape(array, shape)
        reducedIPBKernel = np.sum(array, axis=(0, 1))
        # normalize
        dose = dose / np.max(dose)
        reducedIPBKernel = reducedIPBKernel / np.max(reducedIPBKernel)
        plt.plot(depth, reducedIPBKernel)
        plt.plot(depth, dose)
        plt.xlabel('depth [cm]')
        plt.ylabel('dose (a.u.)')
        plt.legend(['direct simulation', 'Terma to dose'])
        # plt.show()
        outputFile = '/data/qifan/projects/EndtoEnd4/results/' \
            'figures/compare.png'
        plt.title('directly simulated dose v.s. point kernel convolution')
        plt.savefig(outputFile)
        plt.clf()

        plt.plot(depth, Terma)
        plt.xlabel('depth [cm]')
        plt.ylabel('Terma (a.u.)')
        plt.title('Terma v.s. depth')
        outputFile = '/data/qifan/projects/EndtoEnd4/results/' \
            'figures/Terma.png'
        plt.savefig(outputFile)
        plt.clf()


def examinePointOnly():
    filePath = '/data/qifan/projects/EndtoEnd4/' \
        'results/point6MeV1e7Step/array.bin'
    shape = (199, 199, 400)
    array = np.fromfile(filePath, dtype=np.double)
    array = np.reshape(array, shape)
    reducedPointKernel = np.sum(array, axis=(0, 1))
    depth = np.arange(-50, 350) * 0.1  # cm
    plt.plot(depth, reducedPointKernel)
    plt.show()


def countValidEvents():
    """
    This function processes the log files, and count the number 
    of events and number of hits that are at the origin.
    """
    Folder = '/data/qifan/projects/EndtoEnd4/results/point6MeV32'
    template = os.path.join(Folder, 'G4WT*.txt')
    files = glob.glob(template)
    files.sort()

    totalEvents = 0
    totalOrigins = 0
    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()
        context = ''.join(lines)
        localEvents = context.count('Event ID')
        localOrigins = context.count('(0, 0, 0)')
        if (localEvents != localOrigins):
            print('in file {}, the number of total events '
                '!= total origins'.format(file))
        totalEvents += localEvents
        totalOrigins += localOrigins
    # print('total events: {}, total origins: {}'.format(totalEvents, totalOrigins))


if __name__ == '__main__':
    # PartiTrajComp()
    # waterKernelOutputProcessing()
    # examineIPBdose()
    examinePointOnly()
    # countValidEvents()