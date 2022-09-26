import os
import numpy as np
import matplotlib.pyplot as plt


def view_extended_fluence_map():
    FM_dimension = 128
    FM_convolution_radius = 64
    extended_fluence_map_dimension = (FM_dimension + 4 * FM_convolution_radius,
                                      FM_dimension + 4 * FM_convolution_radius)

    extended_fluence_map_path = '/home/qlyu/ShengNAS2/SharedProjectData/' \
                                'QX_beam_orientation/patient1_E2E_output/extended_fluence_map.dat'
    extended_fluence_map = np.fromfile(extended_fluence_map_path, dtype=np.float32)
    extended_fluence_map = np.reshape(extended_fluence_map, extended_fluence_map_dimension)
    print(np.max(extended_fluence_map))
    plt.imshow(extended_fluence_map)
    plt.show()


def view_dose_calculation():
    dose_path = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E_output/dose001.dat'
    dose_shape = (200, 200, 200)
    dose = np.fromfile(dose_path, dtype=np.float32)
    dose = np.reshape(dose, dose_shape)

    azimuth = 0.9424787940998465
    center = np.array((103, 106))

    isocenter_value = dose[103, 106, 81]
    print('isocenter dose value {}'.format(isocenter_value))
    print('maximum dose value {}'.format(np.max(dose)))

    # output_slice = get_cross_section(dose, dose_shape, center, azimuth)
    # plt.imshow(output_slice)
    # plt.show()


def get_cross_section(dose, shape, center, azimuth):
    # dose: 3-dimensional array
    # shape: array<int, 3>
    # center: array<int, 2>

    eps = 1e-7
    norm_direction = (-np.cos(azimuth), -np.sin(azimuth))
    norm_direction = np.array((np.sign(norm_direction[0]) * (abs(norm_direction[0]) + eps),
                      np.sign(norm_direction[1]) * (abs(norm_direction[1]) + eps)))
    k0 = -center[0] / norm_direction[0]
    k1 = -center[1] / norm_direction[1]
    if abs(k0) > abs(k1):
        k0 = k1

    k2 = (shape[0] - center[0]) / norm_direction[0]
    k3 = (shape[1] - center[1]) / norm_direction[1]
    if abs(k2) > abs(k3):
        k2 = k3

    # now, k0 and k2 are two scales
    if k0 < 0 and k2 > 0 :
        endpoint0 = center + k2 * norm_direction
        endpoint1 = center + k0 * norm_direction
    else:
        endpoint0 = center + k0 * norm_direction
        endpoint1 = center + k2 * norm_direction

    vector = endpoint1 - endpoint0
    length = np.sqrt(vector[0] * vector[0] + vector[1] * vector[1])
    steps = int(length)
    step = vector / length

    output_slice = np.zeros((shape[2], steps))
    for i in range(steps):
        point = endpoint0 + i * step
        point00 = (int(np.floor(point[0])), int(np.floor(point[1])))
        weight00 = (point00[0] + 1 - point[0]) * (point00[1] + 1 - point[1])

        point01 = (point00[0], point00[1] + 1)
        weight01 = (point00[0] + 1 - point[0]) * (point[1] - point00[1])

        point10 = (point00[0] + 1, point00[1])
        weight10 = (point[0] - point00[0]) * (point00[1] + 1 - point[1])

        point11 = (point00[0] + 1, point00[1] + 1)
        weight11 = (point[0] - point00[0]) * (point[1] - point00[1])

        if point00[0] >= 0 and point00[1] >= 0 and point00[0] < shape[0] and point00[1] < shape[1]:
            value00 = dose[point00[0], point00[1], :]
        else:
            value00 = 0

        if point01[0] >= 0 and point01[1] >= 0 and point01[0] < shape[0] and point01[1] < shape[1]:
            value01 = dose[point01[0], point01[1], :]
        else:
            value01 = 0

        if point10[0] >= 0 and point10[1] >= 0 and point10[0] < shape[0] and point10[1] < shape[1]:
            value10 = dose[point10[0], point10[1], :]
        else:
            value10 = 0

        if point11[0] >= 0 and point11[1] >= 0 and point11[0] < shape[0] and point11[1] < shape[1]:
            value11 = dose[point11[0], point11[1], :]
        else:
            value11 = 0

        output_slice[:, i] = weight00 * value00 + weight01 * value01 + weight10 * value10 + weight11 * value11
    return output_slice


def view_slices():
    input_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E'
    output_folder = '/data/qifan/projects_qlyu/EndtoEnd3/data'
    items = ['PTV_weight', 'PTV_target', 'OAR_weight', 'OAR_target']
    shape = (200, 200, 197)

    for item in items:
        input_file = os.path.join(input_folder, item + '.dat')
        input_array = np.fromfile(input_file, dtype=np.float32)
        input_array = np.reshape(input_array, shape)

        output_directory = os.path.join(output_folder, item)
        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

        for i in range(shape[2]):
            output_file = os.path.join(output_directory, '{:03d}.png'.format(i+1))
            plt.imsave(output_file, input_array[:, :, i])
        print('{} done!'.format(item))


def view_Dose():
    source0 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary'
    source1 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary_0init'
    source2 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary_smoothness'
    source3 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary_smoothness_eta1e3'
    dest0 = '/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_optimize_stationary'
    dest1 = '/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_optimize_stationary_0init'
    dest2 = '/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_optimize_stationary_smoothness'
    dest3 = '/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_optimize_stationary_smoothness_eta1e3'
    sources = [source0, source1, source2, source3]
    dests = [dest0, dest1, dest2, dest3]
    dose_shape = (200, 200, 200)
    # for i in range(2):
    for i in range(3, 4):
        source = sources[i]
        dest = dests[i]
        if not os.path.isdir(dest):
            os.mkdir(dest)
        source_file = os.path.join(source, 'totalDose.dat')
        dest_folder = os.path.join(dest, 'totalDose')
        if not os.path.isdir(dest_folder):
            os.mkdir(dest_folder)
        dose = np.fromfile(source_file, dtype=np.float32)
        dose = np.reshape(dose, dose_shape)
        max_dose = np.max(dose)
        for j in range(dose_shape[2]):
            dest_file = os.path.join(dest_folder, '{:03d}.png'.format(j+1))
            plt.imsave(dest_file, dose[:, :, j], vmin=0, vmax=max_dose)
        print(i)


def view_loss():
    loss_file_0 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation' \
                '/patient1_optimize_stationary/loss.dat'
    loss_file_1 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation' \
                '/patient1_optimize_stationary_0init/loss.dat'
    loss_file_2 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation' \
                  '/patient1_optimize_stationary_smoothness/DoseLoss.dat'
    loss0 = np.fromfile(loss_file_0, dtype=np.float32)
    loss1 = np.fromfile(loss_file_1, dtype=np.float32)
    loss2 = np.fromfile(loss_file_2, dtype=np.float32)
    assert (len(loss0) == len(loss1))
    concern_range = [40, 100]
    x_axis = np.arange(concern_range[0], concern_range[1])+1
    loss0_range = loss0[concern_range[0]: concern_range[1]]
    loss1_range = loss1[concern_range[0]: concern_range[1]]
    loss2_range = loss2[concern_range[0]: concern_range[1]]
    plt.plot(x_axis, loss0_range)
    plt.plot(x_axis, loss1_range)
    plt.plot(x_axis, loss2_range)
    plt.legend(['non-zero init', 'zero init', 'smoothness'])
    # plt.legend(['non-zero init', 'smoothness'])
    plt.show()

    output_array = np.concatenate([np.expand_dims(loss1_range, axis=1),
                                  np.expand_dims(loss2_range, axis=1)], axis=1)
    print(output_array)


def move_file():
    folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary_0init'
    target = os.path.join(folder, 'fluence_maps')
    num_beams = 20
    for i in range(num_beams):
        source_file = os.path.join(folder, 'fluence_maps{:03d}.dat'.format(i+1))
        target_file = os.path.join(target, '{:03d}.dat'.format(i+1))
        command = 'mv {} {}'.format(source_file, target_file)
        os.system(command)


def view_fluence_map():
    source0 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary'
    source1 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary_0init'
    source3 = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary_smoothness_eta1e3'
    dest0 = '/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_optimize_stationary'
    dest1 = '/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_optimize_stationary_0init'
    dest3 = '/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_optimize_stationary_smoothness_eta1e3'
    sources = [source0, source1, source3]
    dests = [dest0, dest1, dest3]
    num_beams = 20
    fluence_map_shape = (128, 128)
    # for i in range(2):
    for i in range(2, 3):
        source = sources[i]
        dest = dests[i]
        if not os.path.isdir(dest):
            os.mkdir(dest)
        fluence_map_source = os.path.join(source, 'fluence_maps')
        fluence_map_dest = os.path.join(dest, 'fluence_maps')
        if not os.path.isdir(fluence_map_dest):
            os.mkdir(fluence_map_dest)
        for j in range(num_beams):
            fms = os.path.join(fluence_map_source, '{:03d}.dat'.format(j+1))
            fmp = os.path.join(fluence_map_dest, '{:03d}.png'.format(j+1))
            fluence_map = np.fromfile(fms, dtype=np.float32)
            fluence_map = np.reshape(fluence_map, fluence_map_shape)
            plt.imsave(fmp, fluence_map)


def view_water_dose():
    waterDosePath = '/data/qifan/projects_qlyu/EndtoEnd3/data/water_out/waterDose_1.57_0.00.dat'
    shape = (200, 200, 200)
    waterDose = np.fromfile(waterDosePath, dtype=np.float32)
    waterDose = np.reshape(waterDose, shape)

    # central_slice = waterDose[:, :, 100]
    # plt.imshow(central_slice)
    # plt.show()

    center_line = waterDose[:, 100, 100]
    plt.plot(np.arange(shape[0]), center_line)
    plt.show()


if __name__ == '__main__':
    # view_extended_fluence_map()
    # view_dose_calculation()
    # view_slices()
    # view_Dose()
    # view_loss()
    # move_file()
    view_fluence_map()
    # view_loss()
    # view_fluence_map()
    # view_water_dose()