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


def read_annealing_angles():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    iterations = 1000
    num_beams = 20
    num_perturbations = 2
    for i in range(num_patients):
        result_folder = os.path.join(global_folder, 'patient{}_optimize_annealing_1000iters'.format(i+1))
        zenith_file = os.path.join(result_folder, 'zenith.dat')
        azimuth_file = os.path.join(result_folder, 'azimuth.dat')
        taken_file = os.path.join(result_folder, 'taken.dat')
        DoseLoss_file = os.path.join(result_folder, 'DoseLoss.dat')
        SmoothnessLoss_file = os.path.join(result_folder, 'SmoothnessLoss.dat')

        zenith = np.fromfile(zenith_file, dtype=np.float32)
        azimuth = np.fromfile(azimuth_file, dtype=np.float32)
        taken = np.fromfile(taken_file, dtype=bool)
        DoseLoss = np.fromfile(DoseLoss_file, dtype=np.float32)
        SmoothnessLoss = np.fromfile(SmoothnessLoss_file, dtype=np.float32)

        shape0 = (iterations*num_beams, num_perturbations)
        zenith = np.reshape(zenith, shape0)
        azimuth = np.reshape(azimuth, shape0)
        taken = np.int32(taken)
        zenith = zenith[np.arange(shape0[0]), taken]
        azimuth = azimuth[np.arange(shape0[0]), taken]
        shape1 = (iterations, num_beams)
        zenith = np.reshape(zenith, shape1)
        azimuth = np.reshape(azimuth, shape1)

        SmoothnessLoss = np.reshape(SmoothnessLoss, shape1)
        SmoothnessLoss = np.sum(SmoothnessLoss, axis=1)
        totalLoss = DoseLoss + SmoothnessLoss

        DoseLossArgmin = np.argmin(DoseLoss)
        totalLossArgmin = np.argmin(totalLoss)

        # print('for patient {}, minimum dose loss is achieved at {}, minimum total '
        #       'dose is achieved at {}'.format(i+1, DoseLossArgmin, totalLossArgmin))
        # if i == 3:
        #     print(DoseLoss[-20:])

        # use the last beam angle as annealing-selected angles
        zenith_selected = zenith[-1, :]
        azimuth_selected = azimuth[-1, :]
        content = ''
        for j in range(num_beams):
            content = content + '{},{}\n'.format(zenith_selected[j], azimuth_selected[j])

        output_file = os.path.join(global_folder, 'patient{}_E2E'.format(i+1), 'beam_angles_annealing.txt')
        with open(output_file, 'w') as f:
            f.writelines(content)


def convert_zenith_azimuth_to_VarianIEC():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    num_beams = 20

    for i in range(num_patients):
        patient_folder = os.path.join(global_folder, 'patient{}_E2E'.format(i+1))
        annealing_angle_file = os.path.join(patient_folder, 'beam_angles_annealing_correct.txt')
        with open(annealing_angle_file, 'r') as f:
            lines = f.readlines()
        assert len(lines) == num_beams
        zenith_azimuth = np.zeros((num_beams, 2), dtype=np.float32)
        for j in range(num_beams):
            line = lines[j]
            line = line.split(',')
            zenith_azimuth[j, 0] = float(line[0])
            zenith_azimuth[j, 1] = float(line[1])
        # print(zenith_azimuth)
        vector = zenith_azimuth_to_vector(zenith_azimuth)
        gantry, couch = vector_to_VarianIEC(vector)

        # # verification
        # gantry = np.expand_dims(gantry, axis=1)
        # couch = np.expand_dims(couch, axis=1)
        # vector_ = convert(gantry, couch)
        # inner_product = np.sum(vector * vector_, axis=1)
        # print(inner_product)

        content = 'gantryVarianIEC,couchVarianIEC\n'
        for i in range(num_beams):
            content = content + '{},{}\n'.format(gantry[i], couch[i])
        output_file = os.path.join(patient_folder, 'beam_annealing_correct_varianIEC.csv')
        with open(output_file, 'w') as f:
            f.writelines(content)


def convert_zenith_azimuth_to_VarianIEC_additional():
    patient_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E'
    angle_file = os.path.join(patient_folder, 'beam_angles_annealing_correct_additional.txt')
    num_beams = 20
    with open(angle_file, 'r') as f:
        lines = f.readlines()
    assert(len(lines) == num_beams)
    zenith_azimuth = np.zeros((num_beams, 2), dtype=np.float32)
    for i in range(num_beams):
        line = lines[i]
        line = line.split(',')
        zenith_azimuth[i, 0] = float(line[0])
        zenith_azimuth[i, 1] = float(line[1])
    vector = zenith_azimuth_to_vector(zenith_azimuth)
    gantry, couch = vector_to_VarianIEC(vector)
    content = 'gantryVarianIEC,couchVarianIEC\n'
    for i in range(num_beams):
        content = content + '{},{}\n'.format(gantry[i], couch[i])
    output_file = os.path.join(patient_folder, 'beam_annealing_correct_additional_varianIEC.csv')
    with open(output_file, 'w') as f:
        f.writelines(content)


def zenith_azimuth_to_vector(zenith_azimuth):
    zenith = zenith_azimuth[:, 0]
    azimuth = zenith_azimuth[:, 1]
    vector_x = np.sin(zenith) * np.cos(azimuth)
    vector_y = np.sin(zenith) * np.sin(azimuth)
    vector_z = np.cos(zenith)
    vector = np.concatenate((np.expand_dims(vector_x, axis=1),
                             np.expand_dims(vector_y, axis=1),
                             np.expand_dims(vector_z, axis=1)), axis=1)
    return vector


def vector_to_VarianIEC(vector):
    vector_x = vector[:, 0]
    vector_y = vector[:, 1]
    vector_z = vector[:, 2]

    eps = 1e-4
    gantry = np.arccos(vector_x)
    _couch = np.arctan(abs(vector_y) / (abs(vector_z) + eps))
    couch = (vector_y > 0) * (vector_z > 0) * _couch + \
            (vector_y > 0) * (vector_z <=0) * (np.pi - _couch) + \
            (vector_y <= 0) * (vector_z > 0) * (-_couch) + \
            (vector_y <= 0) * (vector_z <= 0) * (np.pi + _couch)
    couch = couch + np.pi / 2

    # rad to degree
    gantry = gantry / np.pi * 180
    couch = couch / np.pi * 180
    return gantry, couch


if __name__ == '__main__':
    # vector_rotation_examination()
    # vector_to_angle_examination()
    # main__()
    # main()
    # read_annealing_angles()
    # convert_zenith_azimuth_to_VarianIEC()
    convert_zenith_azimuth_to_VarianIEC_additional()