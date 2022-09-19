import os
import numpy as np
import matplotlib.pyplot as plt


def view_slices():
    """This function takes as input the folder containing the binary
    files and converts them into images for visualization. The user
    needs to specify the input folder and output directory."""

    # please specify
    input_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E'
    output_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E_output'

    assert os.path.isdir(input_folder), "input folder doesn\'t exist"
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    CT_path = os.path.join(input_folder, 'CT.dat')
    PTV_target_path = os.path.join(input_folder, 'PTV_target.dat')
    PTV_weight_path = os.path.join(input_folder, 'PTV_weight.dat')
    OAR_target_path = os.path.join(input_folder, 'OAR_target.dat')
    OAR_weight_path = os.path.join(input_folder, 'OAR_weight.dat')
    shape_path = os.path.join(input_folder, 'shape.txt')

    with open(shape_path, 'r') as f:
        line = f.readline()
    line = line.split(' ')
    shape = [int(a) for a in line]

    CT = np.fromfile(CT_path, dtype=np.float32)
    PTV_target = np.fromfile(PTV_target_path, dtype=np.float32)
    PTV_weight = np.fromfile(PTV_weight_path, dtype=np.float32)
    OAR_target = np.fromfile(OAR_target_path, dtype=np.float32)
    OAR_weight = np.fromfile(OAR_weight_path, dtype=np.float32)

    CT = np.reshape(CT, shape)
    PTV_target = np.reshape(PTV_target, shape)
    PTV_weight = np.reshape(PTV_weight, shape)
    OAR_target = np.reshape(OAR_target, shape)
    OAR_weight = np.reshape(OAR_weight, shape)

    # # to view
    # print('{} {} {} {}'.format(np.sum(PTV_target), np.sum(PTV_weight), np.sum(OAR_target), np.sum(OAR_weight)))

    CT_output = os.path.join(output_folder, 'CT')
    PTV_weight_output = os.path.join(output_folder, 'PTV_weight')
    OAR_weight_output = os.path.join(output_folder, 'OAR_weight')

    if not os.path.isdir(CT_output):
        os.mkdir(CT_output)
    if not os.path.isdir(PTV_weight_output):
        os.mkdir(PTV_weight_output)
    if not os.path.isdir(OAR_weight_output):
        os.mkdir(OAR_weight_output)

    for i in range(shape[2]):
        output_file_name = os.path.join(CT_output, '{:03d}.png'.format(i+1))
        plt.imsave(output_file_name, CT[:, :, i])
    for i in range(shape[2]):
        output_file_name = os.path.join(PTV_weight_output, '{:03d}.png'.format(i+1))
        plt.imsave(output_file_name, PTV_weight[:, :, i])
    for i in range(shape[2]):
        output_file_name = os.path.join(OAR_weight_output, '{:03d}.png'.format(i+1))
        plt.imsave(output_file_name, OAR_weight[:, :, i])


if __name__ == '__main__':
    view_slices()