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
    plt.imshow(extended_fluence_map)
    plt.show()


def view_dose_calculation():
    dose_path = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E_output/dose006.dat'
    dose_shape = (200, 200, 200)
    dose = np.fromfile(dose_path, dtype=np.float32)
    dose = np.reshape(dose, dose_shape)

    azimuth = -1.361356
    center = np.array((103, 106))
    output_slice = get_cross_section(dose, dose_shape, center, azimuth)
    plt.imshow(output_slice)
    plt.show()


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


if __name__ == '__main__':
    # view_extended_fluence_map()
    view_dose_calculation()