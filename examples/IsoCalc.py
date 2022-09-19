import numpy as np
import os
import matplotlib.pyplot as plt


def write_PTV_binary(PTV_weight, path):
    if os.path.isdir(path):
        os.system('rm -r {}'.format(path))
    os.mkdir(path)
    shape = PTV_weight.shape
    for i in range(shape[2]):
        slice = PTV_weight[:, :, i]
        file_path = os.path.join(path, '{:03}.png'.format(i+1))
        plt.imsave(file_path, slice)


def calculate_isocenter_proprietary():
    """This function calculates the isocenter coordinates as the center of PTV.
    The user should specify voxel size and PTV mask path"""

    target_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    # num_patients = 1
    for i in range(num_patients):
        patient_folder = os.path.join(target_folder, 'patient{}_E2E'.format(i+1))
        PTV_weight_path = os.path.join(patient_folder, 'PTV_weight.dat')
        shape_file = os.path.join(patient_folder, 'shape.txt')
        voxel_size = 2  # unit: mm

        # read shape
        with open(shape_file, 'r') as f:
            line = f.readline()
        line = line.split(' ')
        shape = [int(a) for a in line]

        PTV_weight = np.fromfile(PTV_weight_path, dtype=np.float32)
        PTV_weight = np.reshape(PTV_weight, shape)
        PTV_weight_binary = PTV_weight > 0

        # # write PTV weight
        # PTV_weight_path = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E_output/PTV_weight'
        # write_PTV_binary(PTV_weight_binary, PTV_weight_path)

        x_coordinates = (np.arange(shape[0]) + 0.5) * voxel_size
        x_coordinates = np.expand_dims(x_coordinates, axis=(1, 2))
        y_coordinates = (np.arange(shape[1]) + 0.5) * voxel_size
        y_coordinates = np.expand_dims(y_coordinates, axis=(0, 2))
        z_coordinates = (np.arange(shape[2]) + 0.5) * voxel_size
        z_coordinates = np.expand_dims(z_coordinates, axis=(0, 1))

        num_voxels = np.sum(PTV_weight_binary)
        iso_x = np.sum(PTV_weight_binary * x_coordinates) / num_voxels
        iso_y = np.sum(PTV_weight_binary * y_coordinates) / num_voxels
        iso_z = np.sum(PTV_weight_binary * z_coordinates) / num_voxels
        content = '{},{},{}'.format(iso_x, iso_y, iso_z)
        print(content)

        output_file = os.path.join(patient_folder, 'isocenter.txt')
        with open(output_file, 'w') as f:
            f.writelines(content)


def calculate_isocenter():
    """This function calculates the isocenter coordinates as the center of PTV.
    The user should specify voxel size and PTV mask path"""

    target_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    patient_folder = os.path.join(target_folder, 'patient1_E2E')
    PTV_weight_path = os.path.join(patient_folder, 'PTV_weight.dat')
    shape_file = os.path.join(patient_folder, 'shape.txt')
    voxel_size = 2  # unit: mm

    # read shape
    with open(shape_file, 'r') as f:
        line = f.readline()
    line = line.split(' ')
    shape = [int(a) for a in line]

    PTV_weight = np.fromfile(PTV_weight_path, dtype=np.float32)
    PTV_weight = np.reshape(PTV_weight, shape)
    PTV_weight_binary = PTV_weight > 0

    # # write PTV weight
    # PTV_weight_path = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E_output/PTV_weight'
    # write_PTV_binary(PTV_weight_binary, PTV_weight_path)

    x_coordinates = (np.arange(shape[0]) + 0.5) * voxel_size
    x_coordinates = np.expand_dims(x_coordinates, axis=(1, 2))
    y_coordinates = (np.arange(shape[1]) + 0.5) * voxel_size
    y_coordinates = np.expand_dims(y_coordinates, axis=(0, 2))
    z_coordinates = (np.arange(shape[2]) + 0.5) * voxel_size
    z_coordinates = np.expand_dims(z_coordinates, axis=(0, 1))

    num_voxels = np.sum(PTV_weight_binary)
    iso_x = np.sum(PTV_weight_binary * x_coordinates) / num_voxels
    iso_y = np.sum(PTV_weight_binary * y_coordinates) / num_voxels
    iso_z = np.sum(PTV_weight_binary * z_coordinates) / num_voxels
    content = '{},{},{}'.format(iso_x, iso_y, iso_z)
    print(content)

    # output_file = os.path.join(patient_folder, 'isocenter.txt')
    # with open(output_file, 'w') as f:
    #     f.writelines(content)


def fluence_map_init():
    """This function is to initialize the fluence map for dose calculation.
    The user should specify output folder and fluence map dimension"""

    # please specify
    fluence_map_path = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/fluence_maps'
    fluence_map_dim = (128, 128)
    num_beams = 20

    # this function initializes all fluence maps to all ones
    if os.path.isdir(fluence_map_path):
        os.system('rm -r {}'.format(fluence_map_path))
    os.mkdir(fluence_map_path)

    fluence_map = np.ones(fluence_map_dim, dtype=np.float32)
    for i in range(num_beams):
        output_file_path = os.path.join(fluence_map_path, '{:03d}.dat'.format(i+1))
        fluence_map.tofile(output_file_path)


if __name__ == '__main__':
    calculate_isocenter()
    # fluence_map_init()