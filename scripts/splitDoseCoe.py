import os
import h5py

expFolder = '/media/raid0/qifan/datasets/UCLAPatients/experiment'

def visualize():
    """
    The original hdf5 file is too large, here we try to divide it.
    In summary, the hdf5 data structure is as follows:

    root
    |--calc_specs
    |--filetype
    |--beam
    |  |--metadata
    |  |  |--beam00000
    |  |  |--beam00001
    |  |  | ...
    |  |--data
    |  |  |--beam00000
    |  |  |  |--beamlet_00089
    |  |  |  |  |--coeffs
    |  |  |  |  |--lindex
    |  |  |  |--beamlet_00090
    |  |  |  |  |--coeffs
    |  |  |  |  |--lindex
    |  |  |  | ...
    |  |  |--beam00001
    |  |  |  |--beamlet...
    |  |  |  | ...
    |  |  | ...

    """
    patList = [1]
    for idx in patList:
        patientName = 'patient{}'.format(idx)
        expPatFolder = os.path.join(expFolder, patientName)
        targetFile = os.path.join(expPatFolder, 'Dose_Coefficients.h5')
        file = h5py.File(targetFile, 'r')

        # # view keys
        # # there are three keys, 'beams', 'calc_specs', 'filetype'
        # keys = list(file.keys())
        # print(keys)
        # break

        beams = file['beams']
        calc_specs = file['calc_specs']
        filetype = file['filetype']
        
        # # study filetype
        # # filetype and calc_specs are empty groups
        # # beams has two keys: 'data' and 'metadata'
        # filetypeList = list(filetype)
        # calc_specsList = list(calc_specs)
        # beamsList = list(beams)
        # print(filetypeList)
        # print(calc_specsList)
        # print(beamsList)
        # break
        
        beams_data = beams['data']
        beams_metadata = beams['metadata']
        beams_keys = list(beams_data)
        beams_keys.sort()
        nBeams = len(beams_keys)
        
        for j, beamKey in enumerate(beams_keys):
            beamData = beams_data[beamKey]
            beamMetaData = beams_metadata[beamKey]
            beamDataKeys = list(beamData)
            beamMetaDataKeys = list(beamMetaData)

            # # beamMetaData is empty
            # # beamDataKeys are beamlets
            # print(beamDataKeys, '\n')
            # print(beamMetaDataKeys)
            # break

            for k, beamletKey in enumerate(beamDataKeys):
                beamlet = beamData[beamletKey]
                # # beamlet has two attributes: 'coeffs' and 'lindex'
                # print(list(beamlet))

                coeffs = beamlet['coeffs'][()]
                lindex = beamlet['lindex'][()]
                print(coeffs.size, lindex.size)
                # break
            break
        break


def divide():
    """
    Here we try to divide the original large hdf5 file
    """
    beamsPerFile = 100
    patList = [1]
    for idx in patList:
        patientName = 'patient{}'.format(idx)
        expPatFolder = os.path.join(expFolder, patientName)
        targetFile = os.path.join(expPatFolder, 'Dose_Coefficients.h5')
        file = h5py.File(targetFile, 'r')

        beams = file['beams']
        beams_data = beams['data']
        beams_metadata = beams['metadata']
        beamList = list(beams_data)
        beamList.sort()
        nBeams = len(beamList)

        subFolder = os.path.join(expPatFolder, 'Dose_Coefficients')
        if not os.path.isdir(subFolder):
            os.makedirs(subFolder)
        nFiles = int((nBeams - 1) // beamsPerFile + 1)
        for j in range(nFiles):
        # for j in [5]:
            outFile = os.path.join(subFolder, 'partition{}.h5'.format(j+1))
            if os.path.isfile(outFile):
                os.remove(outFile)
            f = h5py.File(outFile, 'a')
            calc_specs = f.create_group('calc_specs')
            filetype = f.create_group('filetype')
            for k in range(beamsPerFile):
                beamIdx = j * beamsPerFile + k
                if beamIdx >= nBeams:
                    break
                beamName = beamList[beamIdx]
                beamData = beams_data[beamName]
                beamMetaData = beams_metadata[beamName]
                f.create_group(os.path.join('beam', 'metadata', beamName))
                beamletNames = list(beamData)
                beamletNames.sort()
                for beamletName in beamletNames:
                    beamletData = beamData[beamletName]
                    coeffs = beamletData['coeffs'][()]
                    lindex = beamletData['lindex'][()]
                    coeffsPath = os.path.join('beam', 'data', beamName, beamletName, 'coeffs')
                    lindexPath = os.path.join('beam', 'data', beamName, beamletName, 'lindex')
                    f.create_dataset(coeffsPath, data=coeffs)
                    f.create_dataset(lindexPath, data=lindex)
                print(beamName)
            f.close()


if __name__ == '__main__':
    # visualize()
    divide()