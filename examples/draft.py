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


def plot_range(Range, Data, linestyle=None):
    Data_ = Data[Range[0]: Range[1]]
    if linestyle is None:
        plt.plot(np.arange(Range[0], Range[1])+1, Data_)
    else:
        plt.plot(np.arange(Range[0], Range[1])+1, Data_, linestyle)


def compare_stationary_dynamic():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    stationary_folder = 'patient1_optimize_stationary_200iters'
    dynamic_folder = 'patient1_optimize_dynamic'
    stationary_folder = os.path.join(global_folder, stationary_folder)
    dynamic_folder = os.path.join(global_folder, dynamic_folder)

    _DoseLoss = 'DoseLoss.dat'
    _SmoothnessLoss = 'SmoothnessLoss.dat'

    stationaryDoseLossFile = os.path.join(stationary_folder, _DoseLoss)
    stationarySmoothnessLossFile = os.path.join(stationary_folder, _SmoothnessLoss)
    dynamicDoseLossFile = os.path.join(dynamic_folder, _DoseLoss)
    dynamicSmoothnessLossFile = os.path.join(dynamic_folder, _SmoothnessLoss)

    stationaryDoseLoss = np.fromfile(stationaryDoseLossFile, dtype=np.float32)
    stationarySmoothnessLoss = np.fromfile(stationarySmoothnessLossFile, dtype=np.float32)
    dynamicDoseLoss = np.fromfile(dynamicDoseLossFile, dtype=np.float32)
    dynamicSmoothnessLoss = np.fromfile(dynamicSmoothnessLossFile, dtype=np.float32)

    Range = [50, 200]
    plot_range(Range, stationaryDoseLoss)
    plot_range(Range, dynamicDoseLoss)
    plt.legend(['stationaryDoseLoss', 'dynamicDoseLoss'])
    plt.show()
    print('great!')


def calc_relative_angle(angles_org, angles_new):
    """This function calculates the angle difference between the original angle and new angle"""
    assert angles_org.shape == angles_new.shape, "the number of angles are different"
    num_angles = angles_org.shape[0]
    zenith_org = angles_org[:, 0]
    azimuth_org = angles_org[:, 1]
    vector_org = np.zeros((num_angles, 3), dtype=np.float32)
    vector_org[:, 0] = np.sin(zenith_org) * np.cos(azimuth_org)
    vector_org[:, 1] = np.sin(zenith_org) * np.sin(azimuth_org)
    vector_org[:, 2] = np.cos(zenith_org)

    zenith_new = angles_new[:, 0]
    azimuth_new = angles_new[:, 1]
    vector_new = np.zeros((num_angles, 3), dtype=np.float32)
    vector_new[:, 0] = np.sin(zenith_new) * np.cos(azimuth_new)
    vector_new[:, 1] = np.sin(zenith_new) * np.sin(azimuth_new)
    vector_new[:, 2] = np.cos(zenith_new)

    inner_product = vector_org * vector_new
    inner_product = np.sum(inner_product, axis=1)
    angles = np.arccos(inner_product)
    print(inner_product)
    return angles


def angle_comparison():
    angle_org_file = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/beam_angles_E2E.txt'
    zenith_file = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_dynamic/zenith.dat'
    azimuth_file = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_dynamic/azimuth.dat'

    with open(angle_org_file, 'r') as f:
        lines = f.readlines()
    zenith_org = np.zeros((len(lines), 1), dtype=np.float32)
    azimuth_org = np.zeros((len(lines), 1), dtype=np.float32)
    for i, line in enumerate(lines):
        line = line.split(',')
        zenith_org[i, 0] = float(line[0])
        azimuth_org[i, 0] = float(line[1])
    angles_org = np.concatenate((zenith_org, azimuth_org), axis=1)

    iterations = 200
    num_beams = 20
    _zenith = np.fromfile(zenith_file, dtype=np.float32)
    _azimuth = np.fromfile(azimuth_file, dtype=np.float32)
    angles_new = np.zeros((num_beams, 2), dtype=np.float32)
    angles_new[:, 0] = _zenith[(iterations-1)*num_beams: iterations*num_beams]
    angles_new[:, 1] = _azimuth[(iterations-1)*num_beams: iterations*num_beams]

    angles = calc_relative_angle(angles_org, angles_new)
    print(angles)


def view_dynamic_random():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_dynamic_random'
    azimuth_path = os.path.join(global_folder, 'azimuth.dat')
    zenith_path = os.path.join(global_folder, 'zenith.dat')
    angleIndices_path = os.path.join(global_folder, 'angleIndices.dat')

    iterations = 200
    num_beams = 20
    num_perturbations = 5
    azimuth = np.fromfile(azimuth_path, dtype=np.float32)
    zenith = np.fromfile(zenith_path, dtype=np.float32)
    angleIndices = np.fromfile(angleIndices_path, dtype=np.int32)
    assert azimuth.size == iterations * num_beams * num_perturbations
    assert zenith.size == iterations * num_beams * num_perturbations
    assert angleIndices.size == iterations * num_beams

    shape = [iterations * num_beams, num_perturbations]
    zenith = np.reshape(zenith, shape)
    azimuth = np.reshape(azimuth, shape)
    zenith_selected = zenith[np.arange(iterations * num_beams), angleIndices]
    azimuth_selected = azimuth[np.arange(iterations * num_beams), angleIndices]
    zenith_selected = np.reshape(zenith_selected, (iterations, num_beams))
    azimuth_selected = np.reshape(azimuth_selected, (iterations, num_beams))

    for i in range(num_beams):
        plt.plot(np.arange(iterations), zenith_selected[:, i])
    legend = [str(i) for i in range(num_beams)]
    plt.legend(legend)
    plt.title('zenith v.s. iteration')
    plt.show()

    for i in range(num_beams):
        plt.plot(np.arange(iterations), azimuth_selected[:, i])
    plt.legend(legend)
    plt.title('azimuth v.s. iteration')
    plt.show()

    for i in range(num_beams):
        plt.plot(zenith_selected[:, i], azimuth_selected[:, i])
    plt.legend(legend)
    plt.title('azimuth v.s. zenith')
    plt.show()


def compare_dynamic_random():
    dynamic_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_dynamic'
    random_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_dynamic_random'
    iterations = 200
    num_beams = 20
    num_perturbations = 5

    dynamic_zenith_path = os.path.join(dynamic_folder, 'zenith.dat')
    dynamic_azimuth_path = os.path.join(dynamic_folder, 'azimuth.dat')
    dynamic_zenith = np.fromfile(dynamic_zenith_path, dtype=np.float32)
    dynamic_azimuth = np.fromfile(dynamic_azimuth_path, dtype=np.float32)
    dynamic_zenith = np.reshape(dynamic_zenith, (iterations, num_beams))
    dynamic_azimuth = np.reshape(dynamic_azimuth, (iterations, num_beams))

    random_zenith_path = os.path.join(random_folder, 'zenith.dat')
    random_azimuth_path = os.path.join(random_folder, 'azimuth.dat')
    random_angleIndices_path = os.path.join(random_folder, 'angleIndices.dat')
    random_zenith = np.fromfile(random_zenith_path, dtype=np.float32)
    random_azimuth = np.fromfile(random_azimuth_path, dtype=np.float32)
    random_angleIndices = np.fromfile(random_angleIndices_path, dtype=np.int32)
    shape0 = (iterations * num_beams, num_perturbations)
    random_zenith = np.reshape(random_zenith, shape0)
    random_azimuth = np.reshape(random_azimuth, shape0)
    random_zenith_selected = random_zenith[np.arange(shape0[0]), random_angleIndices]
    random_azimuth_selected = random_azimuth[np.arange(shape0[0]), random_angleIndices]
    shape1 = (iterations, num_beams)
    random_zenith_selected = np.reshape(random_zenith_selected, shape1)
    random_azimuth_selected = np.reshape(random_azimuth_selected, shape1)

    # compare
    for i in range(num_beams):
        plt.plot(dynamic_zenith[:, i], dynamic_azimuth[:, i])
    plt.legend([str(i) for i in range(num_beams)])
    plt.title('dynamic')
    plt.show()

    for i in range(num_beams):
        plt.plot(random_zenith_selected[:, i], random_azimuth_selected[:, i])
    plt.legend([str(i) for i in range(num_beams)])
    plt.title('random')
    plt.show()


def generate_uniform_beam_list():
    content = ''
    num_beams = 20
    azimuth_step_size = 2 * np.pi / num_beams
    for i in range(num_beams):
        content = content + '{},{}\n'.format(np.pi/2, azimuth_step_size * i - np.pi)
    dest_file = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/beam_angles_uniform.txt'
    with open(dest_file, 'w') as f:
        f.writelines(content)


def loss_curve_group_meeting():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    group1 = os.path.join(parent_folder, 'patient1_optimize_stationary_200iters')  # stationary optimization, smoothness term
    group2 = os.path.join(parent_folder, 'patient1_optimize_dynamic')  # initialize fluence map with 0
    group3 = os.path.join(parent_folder, 'patient1_optimize_dynamic_random')
    group4 = os.path.join(parent_folder, 'patient1_optimize_dynamic_uniform')
    group_polish = os.path.join(parent_folder, 'patient1_optimize_dynamic_polish')  # initialize with the fluence map of group 1

    groups = [group1, group2, group3, group4]
    iterations = 200
    for i, group in enumerate(groups):
        DoseLoss_path = os.path.join(group, 'DoseLoss.dat')
        DoseLoss = np.fromfile(DoseLoss_path, dtype=np.float32)
        assert DoseLoss.size == iterations
        plot_range([50, 200], DoseLoss)
    legend = ['group{}'.format(i+1) for i in range(len(groups))]
    plt.legend(legend)
    plt.title('loss function')
    plt.xlabel('iterations')
    plt.ylabel('dose loss (a.u.)')
    plt.savefig('/data/qifan/projects_qlyu/EndtoEnd3/data/images/lossComp.png')


def trajectory_group_meeting():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    group1 = os.path.join(parent_folder, 'patient1_optimize_stationary_200iters')  # stationary optimization, smoothness term
    group2 = os.path.join(parent_folder, 'patient1_optimize_dynamic')  # initialize fluence map with 0
    group3 = os.path.join(parent_folder, 'patient1_optimize_dynamic_random')
    group4 = os.path.join(parent_folder, 'patient1_optimize_dynamic_uniform')

    iterations = 200
    num_beams = 20
    num_perturbations = 5

    # for group2
    group2_zenith_path = os.path.join(group2, 'zenith.dat')
    group2_azimuth_path = os.path.join(group2, 'azimuth.dat')
    group2_zenith = np.fromfile(group2_zenith_path, dtype=np.float32)
    group2_azimuth = np.fromfile(group2_azimuth_path, dtype=np.float32)
    group2_zenith = np.reshape(group2_zenith, (iterations, num_beams))
    group2_azimuth = np.reshape(group2_azimuth, (iterations, num_beams))
    for i in range(num_beams):
        plt.plot(group2_zenith[:, i], group2_azimuth[:, i])
    legend = ['beam{}'.format(i+1) for i in range(num_beams)]
    plt.legend(legend)
    plt.xlabel('zenith')
    plt.ylabel('azimuth')
    plt.title('group2')
    # plt.show()
    plt.savefig('/data/qifan/projects_qlyu/EndtoEnd3/data/images/group2_trajectory.png')
    plt.clf()

    # for i in range(num_beams):
    #     plt.plot(np.arange(iterations), group2_zenith[:, i])
    # plt.legend(legend)
    # plt.xlabel('iterations')
    # plt.ylabel('zenith (rad)')
    # plt.savefig()

    group3_zenith_path = os.path.join(group3, 'zenith.dat')
    group3_azimuth_path = os.path.join(group3, 'azimuth.dat')
    group3_angleIndices_path = os.path.join(group3, 'angleIndices.dat')
    group3_zenith = np.fromfile(group3_zenith_path, dtype=np.float32)
    group3_azimuth = np.fromfile(group3_azimuth_path, dtype=np.float32)
    group3_angleIndices = np.fromfile(group3_angleIndices_path, dtype=np.int32)
    shape1 = (iterations * num_beams, num_perturbations)
    group3_zenith = np.reshape(group3_zenith, shape1)
    group3_azimuth = np.reshape(group3_azimuth, shape1)
    group3_zenith_selected = group3_zenith[np.arange(shape1[0]), group3_angleIndices]
    group3_azimuth_selected = group3_azimuth[np.arange(shape1[0]), group3_angleIndices]
    shape2 = (iterations, num_beams)
    group3_zenith_selected = np.reshape(group3_zenith_selected, shape2)
    group3_azimuth_selected = np.reshape(group3_azimuth_selected, shape2)
    for i in range(num_beams):
        plt.plot(group3_zenith_selected[:, i], group3_azimuth_selected[:, i])
    plt.legend(legend)
    plt.xlabel('zenith')
    plt.ylabel('azimuth')
    plt.title('group3')
    # plt.show()
    plt.savefig('/data/qifan/projects_qlyu/EndtoEnd3/data/images/group3_trajectory.png')
    plt.clf()

    group4_zenith_path = os.path.join(group4, 'zenith.dat')
    group4_azimuth_path = os.path.join(group4, 'azimuth.dat')
    group4_zenith = np.fromfile(group4_zenith_path, dtype=np.float32)
    group4_azimuth = np.fromfile(group4_azimuth_path, dtype=np.float32)
    group4_zenith = np.reshape(group4_zenith, (iterations, num_beams))
    group4_azimuth = np.reshape(group4_azimuth, (iterations, num_beams))
    for i in range(num_beams):
        plt.plot(group4_zenith[:, i], group4_azimuth[:, i])
    legend = ['beam{}'.format(i + 1) for i in range(num_beams)]
    plt.legend(legend)
    plt.xlabel('zenith')
    plt.ylabel('azimuth')
    plt.title('group4')
    # plt.show()
    plt.savefig('/data/qifan/projects_qlyu/EndtoEnd3/data/images/group4_trajectory.png')
    plt.clf()


def view_annealing():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_annealing_'
    iterations = 100
    num_beams = 20
    num_perturbations = 2

    perturbationLossPath = os.path.join(global_folder, 'PerturbationLoss.dat')
    perturbationLoss = np.fromfile(perturbationLossPath, dtype=np.float32)
    perturbationLoss = np.reshape(perturbationLoss, (iterations*num_beams, num_perturbations))
    loss_diff = perturbationLoss[:, 0] - perturbationLoss[:, 1]
    loss_diff = np.abs(loss_diff)
    plt.plot(np.arange(iterations*num_beams), loss_diff)
    plt.show()

    interval_interested = loss_diff[-400:]
    diff_mean = np.mean(interval_interested)
    print(diff_mean)

    probabilityPath = os.path.join(global_folder, 'probability.dat')
    probability = np.fromfile(probabilityPath, dtype=np.float32)
    plt.scatter(np.arange(iterations * num_beams), probability)
    plt.show()


def BOO_annealing_compare():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    BOO_folder = os.path.join(parent_folder, 'patient1_optimize_stationary_200iters')
    annealing_folder = os.path.join(parent_folder, 'patient1_optimize_annealing_T1e4_iter200')

    image_output_directory = '/data/qifan/projects_qlyu/EndtoEnd3/data/images'

    iterations = 200
    num_beams = 20
    num_perturbations = 2
    BOO_loss_file = os.path.join(BOO_folder, 'DoseLoss.dat')
    annealing_loss_file = os.path.join(annealing_folder, 'DoseLoss.dat')
    BOO_smoothness_loss_file = os.path.join(BOO_folder, 'SmoothnessLoss.dat')
    annealing_smoothness_loss_file = os.path.join(annealing_folder, 'SmoothnessLoss.dat')
    BOO_loss = np.fromfile(BOO_loss_file, dtype=np.float32)
    annealing_loss = np.fromfile(annealing_loss_file, dtype=np.float32)
    BOO_smoothness_loss = np.fromfile(BOO_smoothness_loss_file, dtype=np.float32)
    annealing_smoothness_loss = np.fromfile(annealing_smoothness_loss_file, dtype=np.float32)
    shape_smoothness = (iterations, num_beams)
    BOO_smoothness_loss = np.reshape(BOO_smoothness_loss, shape_smoothness)
    annealing_smoothness_loss = np.reshape(annealing_smoothness_loss, shape_smoothness)
    BOO_smoothness_loss = np.sum(BOO_smoothness_loss, axis=1)
    annealing_smoothness_loss = np.sum(annealing_smoothness_loss, axis=1)
    Range = [50, 200]
    plot_range(Range, BOO_loss, 'b-')
    plot_range(Range, BOO_smoothness_loss, 'b--')
    plot_range(Range, annealing_loss, 'r-')
    plot_range(Range, annealing_smoothness_loss, 'r--')
    plt.legend(['BOO dose loss', 'BOO smoothness loss', 'annealing dose loss', 'annealing smoothness loss'])
    plt.xlabel('iterations')
    plt.ylabel('loss (a.u.)')
    plt.title('loss comparison')
    # plt.show()
    plt.savefig(os.path.join(image_output_directory, 'annealing_loss_comp.png'))
    plt.clf()

    BOO_total_loss = BOO_loss + BOO_smoothness_loss
    annealing_total_loss = annealing_loss + annealing_smoothness_loss
    plot_range(Range, BOO_total_loss)
    plot_range(Range, annealing_total_loss)
    plt.xlabel('iterations')
    plt.ylabel('loss (a.u.)')
    plt.title('total loss comparison')
    plt.legend(['BOO total loss', 'annealing total loss'])
    # plt.show()
    plt.savefig(os.path.join(image_output_directory, 'annealing_total_loss_comp.png'))
    plt.clf()
    print('the dose loss of BOO algorithm is {}, while that of the '
          'annealing algorithm is {}, the difference is {}'.format(
        BOO_loss[-1], annealing_loss[-1], BOO_loss[-1] - annealing_loss[-1]))
    print('the total loss of BOO algorithm is {}, while that of the '
          'annealing algorithm is {}, the difference is {}'.format(
        BOO_total_loss[-1], annealing_total_loss[-1], BOO_total_loss[-1] - annealing_total_loss[-1]))
    # exit()

    # show beam angle trajectory
    annealing_zenith_file = os.path.join(annealing_folder, 'zenith.dat')
    annealing_azimuth_file = os.path.join(annealing_folder, 'azimuth.dat')
    annealing_taken_file = os.path.join(annealing_folder, 'taken.dat')
    annealing_zenith = np.fromfile(annealing_zenith_file, dtype=np.float32)
    annealing_azimuth = np.fromfile(annealing_azimuth_file, dtype=np.float32)
    annealing_taken = np.fromfile(annealing_taken_file, dtype=bool)
    annealing_taken = np.int32(annealing_taken)
    shape0 = (iterations * num_beams, num_perturbations)
    annealing_zenith = np.reshape(annealing_zenith, shape0)
    annealing_azimuth = np.reshape(annealing_azimuth, shape0)
    annealing_zenith = annealing_zenith[np.arange(shape0[0]), annealing_taken]
    annealing_azimuth = annealing_azimuth[np.arange(shape0[0]), annealing_taken]
    shape1 = (iterations, num_beams)
    annealing_zenith = np.reshape(annealing_zenith, shape1)
    annealing_azimuth = np.reshape(annealing_azimuth, shape1)
    for i in range(num_beams):
        plt.plot(annealing_zenith[:, i], annealing_azimuth[:, i])
    legend = ['beam{}'.format(i+1) for i in range(num_beams)]
    plt.legend(legend)
    plt.xlabel('zenith (rad)')
    plt.ylabel('azimuth (rad)')
    plt.title('azimuth v.s. zenith')
    # plt.show()
    plt.savefig(os.path.join(image_output_directory, 'annealing_azimuth_zenith_plot.png'))
    plt.clf()

    # show angle difference
    for i in range(num_beams):
        plt.plot(np.arange(iterations), annealing_zenith[:, i])
    plt.legend(legend)
    plt.xlabel('iterations')
    plt.ylabel('zenith (rad)')
    plt.title('zenith angle v.s. iterations')
    # plt.show()
    plt.savefig(os.path.join(image_output_directory, 'annealing_zenith_iteration_plot.png'))
    plt.clf()

    for i in range(num_beams):
        plt.plot(np.arange(iterations), annealing_azimuth[:, i])
    plt.legend(legend)
    plt.xlabel('iterations')
    plt.ylabel('azimuth (rad)')
    plt.title('azimuth angle v.s. iterations')
    # plt.show()
    plt.savefig(os.path.join(image_output_directory, 'annealing_azimuth_iteration_plot.png'))
    plt.clf()


def view_annealing_uniform():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    BOO_folder = os.path.join(parent_folder, 'patient1_optimize_stationary_200iters')
    # annealing_uniform_folder = os.path.join(parent_folder, 'patient1_optimize_annealing_uniform')
    annealing_uniform_folder = os.path.join(parent_folder, 'patient1_optimize_annealing_linear')
    dynamic_uniform_folder = os.path.join(parent_folder, 'patient1_optimize_dynamic_uniform')

    iterations = 200
    num_beams = 20
    num_perturbations = 2

    BOO_loss_file = os.path.join(BOO_folder, 'DoseLoss.dat')
    annealing_loss_file = os.path.join(annealing_uniform_folder, 'DoseLoss.dat')
    dynamic_loss_file = os.path.join(dynamic_uniform_folder, 'DoseLoss.dat')
    BOO_loss = np.fromfile(BOO_loss_file, dtype=np.float32)
    annealing_loss = np.fromfile(annealing_loss_file, dtype=np.float32)
    dynamic_loss = np.fromfile(dynamic_loss_file, dtype=np.float32)
    Range = [50, 200]
    plot_range(Range, BOO_loss)
    plot_range(Range, annealing_loss)
    plot_range(Range, dynamic_loss)
    plt.legend(['BOO', 'annealing', 'dynamic'])
    plt.show()

    annealing_zenith = os.path.join(annealing_uniform_folder, 'zenith.dat')
    annealing_azimuth = os.path.join(annealing_uniform_folder, 'azimuth.dat')
    annealing_taken = os.path.join(annealing_uniform_folder, 'taken.dat')
    annealing_zenith = np.fromfile(annealing_zenith, dtype=np.float32)
    annealing_azimuth = np.fromfile(annealing_azimuth, dtype=np.float32)
    annealing_taken = np.fromfile(annealing_taken, dtype=bool)
    annealing_taken = np.int32(annealing_taken)
    shape0 = (iterations * num_beams, num_perturbations)
    annealing_zenith = np.reshape(annealing_zenith, shape0)
    annealing_azimuth = np.reshape(annealing_azimuth, shape0)
    annealing_zenith = annealing_zenith[np.arange(shape0[0]), annealing_taken]
    annealing_azimuth = annealing_azimuth[np.arange(shape0[0]), annealing_taken]
    shape1 = (iterations, num_beams)
    annealing_zenith = np.reshape(annealing_zenith, shape1)
    annealing_azimuth = np.reshape(annealing_azimuth, shape1)

    for i in range(num_beams):
        plt.plot(annealing_zenith[:, i], annealing_azimuth[:, i])
    legend = ['beam{}'.format(i+1) for i in range(num_beams)]
    plt.legend(legend)
    plt.xlabel('zenith')
    plt.ylabel('azimuth')
    plt.title('azimuth v.s. zenith')
    plt.xlim([0, np.pi])
    plt.ylim([-np.pi, np.pi])
    plt.show()

    for i in range(num_beams):
        plt.plot(np.arange(iterations), annealing_zenith[:, i])
    plt.legend(legend)
    plt.xlabel('iterations')
    plt.ylabel('zenith')
    plt.title('zenith v.s. iterations')
    plt.ylim([0, np.pi])
    plt.show()

    for i in range(num_beams):
        plt.plot(np.arange(iterations), annealing_azimuth[:, i])
    plt.legend(legend)
    plt.xlabel('iterations')
    plt.ylabel('azimuth')
    plt.title('azimuth v.s. iterations')
    plt.ylim([-np.pi, np.pi])
    plt.show()


def view_annealing_scheduling():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    stationary_folder = os.path.join(parent_folder, 'patient1_optimize_stationary_200iters')
    linear_folder = os.path.join(parent_folder, 'patient1_optimize_annealing_linear')

    iterations = 200
    num_beams = 20
    num_perturbations = 2
    stationaryDoseLossFile = os.path.join(stationary_folder, 'DoseLoss.dat')
    stationarySmoothnessFile = os.path.join(stationary_folder, 'SmoothnessLoss.dat')
    linearDoseLossFile = os.path.join(linear_folder, 'DoseLoss.dat')
    linearSmoothnessFile = os.path.join(linear_folder, 'SmoothnessLoss.dat')

    stationaryDoseLoss = np.fromfile(stationaryDoseLossFile, dtype=np.float32)
    stationarySmoothness = np.fromfile(stationarySmoothnessFile, dtype=np.float32)
    linearDoseLoss = np.fromfile(linearDoseLossFile, dtype=np.float32)
    linearSmoothness = np.fromfile(linearSmoothnessFile, dtype=np.float32)

    stationarySmoothness = np.reshape(stationarySmoothness, (iterations,num_beams))
    linearSmoothness = np.reshape(linearSmoothness, (iterations, num_beams))
    stationarySmoothness = np.sum(stationarySmoothness, axis=1)
    linearSmoothness = np.sum(linearSmoothness, axis=1)

    Range = [50, 200]
    plot_range(Range, stationaryDoseLoss)
    plot_range(Range, linearDoseLoss)
    plt.legend(['BOO', 'annealing'])
    plt.show()

    annealing_zenith = os.path.join(linear_folder, 'zenith.dat')
    annealing_azimuth = os.path.join(linear_folder, 'azimuth.dat')
    annealing_taken = os.path.join(linear_folder, 'taken.dat')
    annealing_zenith = np.fromfile(annealing_zenith, dtype=np.float32)
    annealing_azimuth = np.fromfile(annealing_azimuth, dtype=np.float32)
    annealing_taken = np.fromfile(annealing_taken, dtype=bool)
    annealing_taken = np.int32(annealing_taken)
    shape0 = (iterations * num_beams, num_perturbations)
    annealing_zenith = np.reshape(annealing_zenith, shape0)
    annealing_azimuth = np.reshape(annealing_azimuth, shape0)
    annealing_zenith = annealing_zenith[np.arange(shape0[0]), annealing_taken]
    annealing_azimuth = annealing_azimuth[np.arange(shape0[0]), annealing_taken]
    shape1 = (iterations, num_beams)
    annealing_zenith = np.reshape(annealing_zenith, shape1)
    annealing_azimuth = np.reshape(annealing_azimuth, shape1)

    for i in range(num_beams):
        plt.plot(annealing_zenith[:, i], annealing_azimuth[:, i])
    legend = ['beam{}'.format(i+1) for i in range(num_beams)]
    plt.legend(legend)
    plt.xlabel('zenith')
    plt.ylabel('azimuth')
    plt.title('azimuth v.s. zenith')
    plt.xlim([0, np.pi])
    plt.ylim([-np.pi, np.pi])
    plt.show()

    for i in range(num_beams):
        plt.plot(np.arange(iterations), annealing_zenith[:, i])
    plt.legend(legend)
    plt.xlabel('iterations')
    plt.ylabel('zenith')
    plt.title('zenith v.s. iterations')
    plt.ylim([0, np.pi])
    plt.show()

    for i in range(num_beams):
        plt.plot(np.arange(iterations), annealing_azimuth[:, i])
    plt.legend(legend)
    plt.xlabel('iterations')
    plt.ylabel('azimuth')
    plt.title('azimuth v.s. iterations')
    plt.ylim([-np.pi, np.pi])
    plt.show()


if __name__ == '__main__':
    # view_extended_fluence_map()
    # view_dose_calculation()
    # view_slices()
    # view_Dose()
    # view_loss()
    # move_file()
    # view_fluence_map()
    # view_loss()
    # view_fluence_map()
    # view_water_dose()
    # compare_stationary_dynamic()
    # angle_comparison()
    # view_dynamic_random()
    # compare_dynamic_random()
    # generate_uniform_beam_list()
    # loss_curve_group_meeting()
    # trajectory_group_meeting()
    # view_annealing()
    # BOO_annealing_compare()
    # view_annealing_uniform()
    view_annealing_scheduling()