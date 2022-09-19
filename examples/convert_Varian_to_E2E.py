import os
import numpy as np


def convert(gantry, couch):
    """This function converts the angles in Varian IEC definition (gantry, couch)
    to the unit vector after rotation. It takes as input two column arrays and
    output an array of unit vector. coordinate framework definition is as follows:
        x: from anterior to posterior
        y: from right to left
        z: from superior to inferior
    Input is in degree, and output is in radian
    """

    gantry_rad = gantry * np.pi / 180
    couch_rad = couch * np.pi / 180

    # initial angle is in (1, 0, 0)
    vector_after_gantry_rotation = np.concatenate((
        np.cos(gantry_rad), -np.sin(gantry_rad), np.zeros_like(gantry_rad)), axis=1)

    vector_after_couch_rotation_0 = vector_after_gantry_rotation[:, 0]
    vector_after_couch_rotation_1 = vector_after_gantry_rotation[:, 1] * np.squeeze(np.cos(couch_rad))
    vector_after_couch_rotation_2 = - vector_after_gantry_rotation[:, 1] * np.squeeze(np.sin(couch_rad))
    vector_after_couch_rotation = np.concatenate((
        np.expand_dims(vector_after_couch_rotation_0, axis=1),
        np.expand_dims(vector_after_couch_rotation_1, axis=1),
        np.expand_dims(vector_after_couch_rotation_2, axis=1)), axis=1)

    return vector_after_couch_rotation


def vector_to_angle(vector):
    """This function converts the unit vector to (zenith, azimuth)"""
    eps = 1e-8
    vector_norm = np.sqrt(vector[:, 0]**2 + vector[:, 1]**2 + vector[:, 2]**2)
    zenith = np.arccos(vector[:, 2] / (vector_norm + eps))
    xy_norm = np.sqrt(vector[:, 0]**2 + vector[:, 1]**2)
    azimuth = np.arccos(vector[:, 0] / (xy_norm + eps)) * np.sign(vector[:, 1])
    return zenith, azimuth


def angle_to_vector(zenith, azimuth):
    vector_after_zenith_0 = np.sin(zenith)
    # vector_after_zenith_1 = np.zeros_like(vector_after_zenith_0)
    vector_after_zenith_2 = np.cos(zenith)

    vector_after_azimuth_0 = vector_after_zenith_0 * np.cos(azimuth)
    vector_after_azimuth_1 = vector_after_zenith_0 * np.sin(azimuth)
    vector_after_azimuth_2 = vector_after_zenith_2

    vector_after_azimuth = np.concatenate((
        np.expand_dims(vector_after_azimuth_0, axis=1),
        np.expand_dims(vector_after_azimuth_1, axis=1),
        np.expand_dims(vector_after_azimuth_2, axis=1)), axis=1)
    return vector_after_azimuth


def vector_rotation_examination():
    angles = np.array([[0, 0], [90, 0], [90, 90], [0, 90], [30, 30]])
    gantry = angles[:, 0]
    couch = angles[:, 1]
    gantry = np.expand_dims(gantry, axis=1)
    couch = np.expand_dims(couch, axis=1)
    vector = convert(gantry, couch)

    print(angles)
    print(vector)


def vector_to_angle_examination():
    # firstly, we generate a random set of (zenith, azimuth) angle pairs
    num_pairs = 16
    zenith_org = np.random.rand(num_pairs) * np.pi
    azimuth_org = (np.random.rand(num_pairs) - 0.5) * 2 * np.pi
    vector_org = angle_to_vector(zenith_org, azimuth_org)

    zenith_new, azimuth_new = vector_to_angle(vector_org)
    print(zenith_org)
    print(zenith_new)
    print(azimuth_org)
    print(azimuth_new)


def main__():
    """This function takes as input the path to the file containing Varian IEC
    angles (gantry,couch), and outputs the file containing (zenith,azimuth).
    The user should specify the input path and the output path"""

    num_patients = 6
    target_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    for idx in range(num_patients):
        patient_folder = os.path.join(target_folder, 'patient{}_E2E'.format(idx+1))
        input_file = os.path.join(patient_folder, 'beam_angles_VarianIEC.csv')
        output_file = os.path.join(patient_folder, 'beam_angles_E2E.txt')
        with open(input_file, 'r') as f:
            lines = f.readlines()
        lines = lines[1:]
        num_beams = len(lines)
        gantry = np.zeros((num_beams, 1))
        couch = np.zeros((num_beams, 1))
        for j, line in enumerate(lines):
            line = line.split(',')
            gantry[j, 0] = float(line[0])
            couch[j, 0] = float(line[1])
        vector_after_couch_rotation = convert(gantry, couch)
        zenith, azimuth = vector_to_angle(vector_after_couch_rotation)

        # convert data to text
        content = ''
        for j in range(zenith.shape[0]):
            content = content + '{},{}\n'.format(zenith[j], azimuth[j])
        with open(output_file, 'w') as f:
            f.writelines(content)


def main():
    """This function takes as input the path to the file containing Varian IEC
    angles (gantry,couch), and outputs the file containing (zenith,azimuth).
    The user should specify the input path and the output path"""

    # please specify
    patient_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E'
    input_file = os.path.join(patient_folder, 'beam_angles_VarianIEC.csv')
    output_file = os.path.join(patient_folder, 'beam_angles_E2E.txt')

    with open(input_file, 'r') as f:
        lines = f.readlines()
    lines = lines[1:]
    num_beams = len(lines)
    gantry = np.zeros((num_beams, 1))
    couch = np.zeros((num_beams, 1))
    for j, line in enumerate(lines):
        line = line.split(',')
        gantry[j, 0] = float(line[0])
        couch[j, 0] = float(line[1])
    vector_after_couch_rotation = convert(gantry, couch)
    zenith, azimuth = vector_to_angle(vector_after_couch_rotation)

    # convert data to text
    content = ''
    for j in range(zenith.shape[0]):
        content = content + '{},{}\n'.format(zenith[j], azimuth[j])
    with open(output_file, 'w') as f:
        f.writelines(content)


if __name__ == '__main__':
    # vector_rotation_examination()
    # vector_to_angle_examination()
    # main__()
    main()