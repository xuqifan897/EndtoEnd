"""
This python script compares the number of particle 
entries and number of tracks. We assume they are equal
"""

import os


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
    outputFolder = './output'
    outputFile = os.path.join(outputFolder, 'myOutput.txt')
    with open(outputFile, 'r') as f:
        lines = f.readlines()
    logStartLine = 947
    logEndLine = 3373
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


if __name__ == '__main__':
    # PartiTrajComp()
    waterKernelOutputProcessing()