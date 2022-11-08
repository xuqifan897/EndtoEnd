import os
import numpy as np
import matplotlib.pyplot as plt
import mat73
import cv2


def draw_annealing_loss():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/'
    num_patients = 6
    iterations = 1000
    num_beams = 20
    num_perturbations = 2
    startPoint = 200

    fig, axes = plt.subplots(2, 3)
    coords = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    for i in range(num_patients):
        if i == 0:
            continue
        optimize_annealing_folder = os.path.join(global_folder, 'patient{}_annealing_correct'.format(i+1))
        doseLossFile = os.path.join(optimize_annealing_folder, 'DoseLoss.dat')
        SmoothnessLossFile = os.path.join(optimize_annealing_folder, 'SmoothnessLoss.dat')
        zenithFile = os.path.join(optimize_annealing_folder, 'zenith.dat')
        azimuthFile = os.path.join(optimize_annealing_folder, 'azimuth.dat')
        takenFile = os.path.join(optimize_annealing_folder, 'taken.dat')

        doseLoss = np.fromfile(doseLossFile, dtype=np.float32)
        SmoothnessLoss = np.fromfile(SmoothnessLossFile, dtype=np.float32)
        zenith = np.fromfile(zenithFile, dtype=np.float32)
        azimuth = np.fromfile(azimuthFile, dtype=np.float32)
        taken = np.fromfile(takenFile, bool)
        taken = np.int32(taken)

        shape0 = (iterations, num_beams)
        SmoothnessLoss = np.reshape(SmoothnessLoss, shape0)
        SmoothnessLoss = np.sum(SmoothnessLoss, axis=1)

        shape0 = (iterations * num_beams, num_perturbations)
        zenith = np.reshape(zenith, shape0)
        azimuth = np.reshape(azimuth, shape0)
        zenith = zenith[np.arange(shape0[0]), taken]
        azimuth = azimuth[np.arange(shape0[0]), taken]

        # plt.plot(np.arange(startPoint, iterations), doseLoss[startPoint: iterations])
        # plt.xlabel('iterations')
        # plt.ylabel('dose loss (a.u.)')
        # plt.ticklabel_format(axis='y', scilimits=[0, 0])
        # plt.title('patient {}'.format(i+1))
        # plt.show()

        coord = coords[i]
        axes[coord[0], coord[1]].plot(np.arange(startPoint, iterations), doseLoss[startPoint: iterations])
        axes[coord[0], coord[1]].set_title('patient {}'.format(i+1), fontsize=20)
        axes[coord[0], coord[1]].ticklabel_format(axis='y', scilimits=[0, 0])
        axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
        axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
        axes[coord[0], coord[1]].yaxis.offsetText.set_fontsize(14)

    # additional experiment on patient 1
    optimize_annealing_folder = os.path.join(global_folder, 'patient11_annealing_correct')
    doseLossFile = os.path.join(optimize_annealing_folder, 'DoseLoss.dat')
    doseLoss = np.fromfile(doseLossFile, dtype=np.float32)
    # plt.plot(np.arange(startPoint, iterations), doseLoss[startPoint: iterations])
    # plt.xlabel('iterations')
    # plt.ylabel('dose loss (a.u.)')
    # plt.ticklabel_format(axis='y', scilimits=[0, 0])
    # plt.title('patient {}'.format(1))
    # plt.show()
    axes[0, 0].plot(np.arange(startPoint, iterations), doseLoss[startPoint: iterations])
    axes[0, 0].set_title('patient {}'.format(1), fontsize=20)
    axes[0, 0].ticklabel_format(axis='y', scilimits=[0, 0])
    axes[0, 0].tick_params(axis='both', which='major', labelsize=14)
    axes[0, 0].tick_params(axis='both', which='major', labelsize=14)
    axes[0, 0].yaxis.offsetText.set_fontsize(14)

    for ax in axes.flat:
        ax.set_xlabel('iterations', fontsize=16)
        ax.set_ylabel('dose loss (a.u.)', fontsize=16)
    # for ax in axes.flat:
    #     ax.label_outer()
    fig.set_figheight(10)
    fig.set_figwidth(20)
    fig.tight_layout()
    # fig.show()
    fig.savefig('./annealingLoss.pdf')


def draw_compare_E2E():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    iterations = 200
    num_beams = 20
    startPoint = 50
    flag = 'dose'
    # flag = 'total'
    roof_scale = 4

    fig, axes = plt.subplots(2, 3)
    coords = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    for i in range(num_patients):
        if i == 0:
            continue
        annealing_init_folder = os.path.join(global_folder, 'patient{}_annealing_correct_init'.format(i+1))
        BOO_init_folder = os.path.join(global_folder, 'patient{}_BOO_correct'.format(i+1))

        annealingDoseLossFile = os.path.join(annealing_init_folder, 'DoseLoss.dat')
        annealingSmoothnessLossFile = os.path.join(annealing_init_folder, 'SmoothnessLoss.dat')
        BOODoseLossFile = os.path.join(BOO_init_folder, 'DoseLoss.dat')
        BOOSmoothnessLossFile = os.path.join(BOO_init_folder, 'SmoothnessLoss.dat')

        annealingDoseLoss = np.fromfile(annealingDoseLossFile, dtype=np.float32)
        annealingSmoothnessLoss = np.fromfile(annealingSmoothnessLossFile, dtype=np.float32)
        BOODoseLoss = np.fromfile(BOODoseLossFile, dtype=np.float32)
        BOOSmoothnessLoss = np.fromfile(BOOSmoothnessLossFile, dtype=np.float32)

        shape = (iterations, num_beams)
        annealingSmoothnessLoss = np.reshape(annealingSmoothnessLoss, shape)
        BOOSmoothnessLoss = np.reshape(BOOSmoothnessLoss, shape)
        annealingSmoothnessLoss = np.sum(annealingSmoothnessLoss, axis=1)
        BOOSmoothnessLoss = np.sum(BOOSmoothnessLoss, axis=1)

        if flag == 'dose':
            interval_annealing = annealingDoseLoss[startPoint:]
            interval_BOO = BOODoseLoss[startPoint:]
            roof = roof_scale * max(np.max(interval_annealing), np.max(interval_BOO))
            # plt.plot(np.arange(iterations), BOODoseLoss)
            # plt.plot(np.arange(iterations), annealingDoseLoss)
            # plt.ylim([0, roof])
            # plt.legend(['BOO', 'annealing'])
            # plt.xlabel('iterations')
            # plt.ylabel('dose loss (a.u.)')
            # plt.title('patient {}'.format(i+1))
            # plt.show()

            coord = coords[i]
            axes[coord[0], coord[1]].plot(np.arange(iterations), BOODoseLoss)
            axes[coord[0], coord[1]].plot(np.arange(iterations), annealingDoseLoss)
            axes[coord[0], coord[1]].set_ylim([0, roof])
            axes[coord[0], coord[1]].legend(['baseline', 'annealing'], fontsize=16)
            axes[coord[0], coord[1]].set_title('patient {}'.format(i+1), fontsize=20)
            axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
            axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
            axes[coord[0], coord[1]].yaxis.offsetText.set_fontsize(14)
        else:
            annealingLoss = annealingDoseLoss + annealingSmoothnessLoss
            BOOLoss = BOODoseLoss + BOOSmoothnessLoss
            interval_annealing = annealingLoss[startPoint:]
            interval_BOO = BOOLoss[startPoint:]
            roof = roof_scale * max(np.max(interval_annealing), np.max(interval_BOO))
            plt.plot(np.arange(iterations), BOOLoss)
            plt.plot(np.arange(iterations), annealingLoss)
            plt.ylim([0, roof])
            plt.legend(['BOO', 'annealing'])
            plt.xlabel('iterations')
            plt.ylabel('total loss (a.u.)')
            plt.title('patient {}'.format(i+1))
            plt.show()

    # additional experiment on patient 1
    annealing_init_folder = os.path.join(global_folder, 'patient1_annealing_correct_additional_init')
    BOO_init_folder = os.path.join(global_folder, 'patient1_BOO_correct')

    annealingDoseLossFile = os.path.join(annealing_init_folder, 'DoseLoss.dat')
    annealingSmoothnessLossFile = os.path.join(annealing_init_folder, 'SmoothnessLoss.dat')
    BOODoseLossFile = os.path.join(BOO_init_folder, 'DoseLoss.dat')
    BOOSmoothnessLossFile = os.path.join(BOO_init_folder, 'SmoothnessLoss.dat')

    annealingDoseLoss = np.fromfile(annealingDoseLossFile, dtype=np.float32)
    annealingSmoothnessLoss = np.fromfile(annealingSmoothnessLossFile, dtype=np.float32)
    BOODoseLoss = np.fromfile(BOODoseLossFile, dtype=np.float32)
    BOOSmoothnessLoss = np.fromfile(BOOSmoothnessLossFile, dtype=np.float32)

    shape = (iterations, num_beams)
    annealingSmoothnessLoss = np.reshape(annealingSmoothnessLoss, shape)
    BOOSmoothnessLoss = np.reshape(BOOSmoothnessLoss, shape)
    annealingSmoothnessLoss = np.sum(annealingSmoothnessLoss, axis=1)
    BOOSmoothnessLoss = np.sum(BOOSmoothnessLoss, axis=1)

    if flag == 'dose':
        interval_annealing = annealingDoseLoss[startPoint:]
        interval_BOO = BOODoseLoss[startPoint:]
        roof = roof_scale * max(np.max(interval_annealing), np.max(interval_BOO))
        # plt.plot(np.arange(iterations), BOODoseLoss)
        # plt.plot(np.arange(iterations), annealingDoseLoss)
        # plt.ylim([0, roof])
        # plt.legend(['BOO', 'annealing'])
        # plt.xlabel('iterations')
        # plt.ylabel('dose loss (a.u.)')
        # plt.title('patient 1')
        # plt.show()

        coord = [0, 0]
        axes[coord[0], coord[1]].plot(np.arange(iterations), BOODoseLoss)
        axes[coord[0], coord[1]].plot(np.arange(iterations), annealingDoseLoss)
        axes[coord[0], coord[1]].set_ylim([0, roof])
        axes[coord[0], coord[1]].legend(['baseline', 'annealing'], fontsize=16)
        axes[coord[0], coord[1]].set_title('patient {}'.format(0), fontsize=20)
        axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
        axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
        axes[coord[0], coord[1]].yaxis.offsetText.set_fontsize(14)
    else:
        annealingLoss = annealingDoseLoss + annealingSmoothnessLoss
        BOOLoss = BOODoseLoss + BOOSmoothnessLoss
        interval_annealing = annealingLoss[startPoint:]
        interval_BOO = BOOLoss[startPoint:]
        roof = roof_scale * max(np.max(interval_annealing), np.max(interval_BOO))
        plt.plot(np.arange(iterations), BOOLoss)
        plt.plot(np.arange(iterations), annealingLoss)
        plt.ylim([0, roof])
        plt.legend(['BOO', 'annealing'])
        plt.xlabel('iterations')
        plt.ylabel('total loss (a.u.)')
        plt.title('patient 1')
        plt.show()

    for ax in axes.flat:
        ax.set_xlabel('iterations', fontsize=16)
        ax.set_ylabel('dose loss (a.u.)', fontsize=16)

    fig.set_figheight(10)
    fig.set_figwidth(20)
    fig.tight_layout()
    # fig.show()
    fig.savefig('./polishingE2E.pdf')


def draw_compare_BOO():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    iterations = 2000
    start = 20
    stop = 100
    flag = 'dose'
    # flag = 'total'

    fig, axes = plt.subplots(2, 3)
    coords = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    for i in range(num_patients):
        patient_folder = os.path.join(parent_folder, 'patient{}_compare'.format(i+1))
        BOO_file = os.path.join(patient_folder, 'optimize_BOO', 'polish_result.mat')
        annealing_file = os.path.join(patient_folder, 'optimize_annealing', 'polish_result.mat')
        if i == 0:
            annealing_file = os.path.join(patient_folder, 'optimize_annealing_additional', 'polish_result.mat')
        BOO_mat = mat73.loadmat(BOO_file)
        annealing_mat = mat73.loadmat(annealing_file)

        # plot costsDF, i.e., dose loss
        if flag == 'dose':
            BOO_costsDF = BOO_mat['result']['costsDF_polish']
            annealing_costsDF = annealing_mat['result']['costsDF_polish']
            # plt.plot(np.arange(start, stop), BOO_costsDF[start:stop])
            # plt.plot(np.arange(start, stop), annealing_costsDF[start:stop])
            # plt.legend(['BOO', 'annealing'])
            # plt.xlabel('iterations')
            # plt.ylabel('dose loss (a.u.)')
            # plt.ticklabel_format(axis='y', scilimits=[0, 0])
            # plt.title('patient {}'.format(i+1))
            # plt.show()

            coord = coords[i]
            axes[coord[0], coord[1]].plot(np.arange(start, stop), BOO_costsDF[start:stop])
            axes[coord[0], coord[1]].plot(np.arange(start, stop), annealing_costsDF[start:stop])
            axes[coord[0], coord[1]].legend(['baseline', 'annealing'], fontsize=16)
            axes[coord[0], coord[1]].set_title('patient {}'.format(0), fontsize=20)
            axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
            axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
            axes[coord[0], coord[1]].yaxis.offsetText.set_fontsize(14)
            axes[coord[0], coord[1]].ticklabel_format(axis='y', scilimits=[0, 0])
            axes[coord[0], coord[1]].set_xlabel('iterations', fontsize=16)
            axes[coord[0], coord[1]].set_ylabel('does loss (a.u.)', fontsize=16)

        # plot costs, i.e., total loss
        else:
            BOO_costs = BOO_mat['result']['costs_polish']
            annealing_costs = annealing_mat['result']['costs_polish']
            plt.plot(np.arange(start, stop), BOO_costs[start:stop])
            plt.plot(np.arange(start, stop), annealing_costs[start:stop])
            plt.legend(['BOO', 'annealing'])
            plt.xlabel('iterations')
            plt.ylabel('total loss (a.u.)')
            plt.ticklabel_format(axis='y', scilimits=[0, 0])
            plt.title('patient {}'.format(i+1))
            plt.show()

    fig.set_figheight(10)
    fig.set_figwidth(20)
    fig.tight_layout()
    # fig.show()
    fig.savefig('./polishingBOO.pdf')


def draw_DVH():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    shapes_org = [(200, 200, 197), (200, 200, 138), (171, 171, 112), (170, 170, 113), (200, 200, 152), (260, 260, 128)]
    num_patients = 6
    append_norm = 8
    prescription_dose = 20
    shapes_append = []
    for i in range(num_patients):
        shape_org = shapes_org[i]
        shape_append = (shape_org[0], shape_org[1], int(append_norm * np.ceil(shape_org[2] / append_norm)))
        shapes_append.append(shape_append)

    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    exclude_list = ['BODY.dat', 'Skin.dat', 'RingStructure.dat']
    fig, axes = plt.subplots(2, 3)
    coords = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    for i in range(num_patients):
        annealing_init_folder = os.path.join(parent_folder, 'patient{}_annealing_correct_init'.format(i+1))
        if i == 0:
            annealing_init_folder = os.path.join(parent_folder, 'patient1_annealing_correct_additional_init')
        totalDoseAnnealingFile = os.path.join(annealing_init_folder, 'totalDose.dat')
        totalDoseAnnealing = np.fromfile(totalDoseAnnealingFile, dtype=np.float32)
        shape_org = shapes_org[i]
        shape_append = shapes_append[i]
        totalDoseAnnealing = np.reshape(totalDoseAnnealing, shape_append)

        BOO_init_folder = os.path.join(parent_folder, 'patient{}_BOO_correct'.format(i+1))
        totalDoseBOOFile = os.path.join(BOO_init_folder, 'totalDose.dat')
        totalDoseBOO = np.fromfile(totalDoseBOOFile, dtype=np.float32)
        totalDoseBOO = np.reshape(totalDoseBOO, shape_append)

        mask_folder = os.path.join(parent_folder, 'patient{}_compare'.format(i+1), 'masks')
        mask_files = os.listdir(mask_folder)
        PTV_file = os.path.join(mask_folder, 'PTV.dat')
        PTV_mask = np.fromfile(PTV_file, dtype=bool)
        PTV_mask = mask_append(PTV_mask, shape_org, shape_append)

        # normalize
        normalized_dose_annealing = dose_normalize(totalDoseAnnealing, PTV_mask)
        normalized_dose_BOO = dose_normalize(totalDoseBOO, PTV_mask)
        totalDoseAnnealing *= prescription_dose / normalized_dose_annealing
        totalDoseBOO *= prescription_dose / normalized_dose_BOO

        # draw DVH
        coord = coords[i]
        for j, mask_file in enumerate(mask_files):
            if mask_file in exclude_list:
                continue
            mask_file = os.path.join(mask_folder, mask_file)
            mask_array = np.fromfile(mask_file, dtype=bool)
            mask_array = mask_append(mask_array, shape_org, shape_append)
            draw_line(totalDoseAnnealing, mask_array, axes[coord[0], coord[1]], color_list[j], '-')

        for j, mask_file in enumerate(mask_files):
            if mask_file in exclude_list:
                continue
            mask_file = os.path.join(mask_folder, mask_file)
            mask_array = np.fromfile(mask_file, dtype=bool)
            mask_array = mask_append(mask_array, shape_org, shape_append)
            draw_line(totalDoseBOO, mask_array, axes[coord[0], coord[1]], color_list[j], '--')

        legend = []
        for mask_file in mask_files:
            if mask_file in exclude_list:
                continue
            structure = mask_file.split('.')[0]
            legend.append(structure)
        axes[coord[0], coord[1]].legend(legend, loc='upper right', fontsize=16)
        axes[coord[0], coord[1]].set_xlabel('dose (Gy)', fontsize=16)
        axes[coord[0], coord[1]].set_ylabel('fractional volume', fontsize=16)
        axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
        axes[coord[0], coord[1]].tick_params(axis='both', which='major', labelsize=14)
        axes[coord[0], coord[1]].set_title('patient {}'.format(i+1), fontsize=20)

    fig.set_figheight(10)
    fig.set_figwidth(20)
    fig.tight_layout()
    # fig.show()
    fig.savefig('./DVH.pdf')


def draw_line(totalDose, mask_array, ax, color, linestyle):
    mask_dose = totalDose[mask_array]
    mask_dose.sort()
    mask_dose = np.insert(mask_dose, 0, 0.)
    percentile_array = np.arange(mask_dose.size, 0, -1) / (mask_dose.size)
    ax.plot(mask_dose, percentile_array, color=color, linestyle=linestyle)


def mask_append(mask_input, shape_org, shape_append):
    mask_output = np.zeros(shape_append, dtype=bool)
    mask_input = np.reshape(mask_input, shape_org)
    mask_output[:, :, :shape_org[2]] = mask_input
    return mask_output


def dose_normalize(totalDose, PTV_mask):
    PTVDose = totalDose[PTV_mask]
    thresh = np.percentile(PTVDose, 5)
    return thresh


def draw_masks():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    data_folder = '/data/qifan/projects_qlyu/EndtoEnd3/data'
    shapes_org = [(200, 200, 197), (200, 200, 138), (171, 171, 112), (170, 170, 113), (200, 200, 152), (260, 260, 128)]
    num_patients = 6
    exclude_list = ['PTV.dat', 'BODY.dat', 'Skin.dat', 'RingStructure.dat']

    for i in range(num_patients):
        shape_org = shapes_org[i]
        OAR_mask = np.zeros(shape_org, dtype=bool)
        mask_folder = os.path.join(parent_folder, 'patient{}_compare'.format(i+1), 'masks')
        mask_files = os.listdir(mask_folder)
        for mask_file in mask_files:
            if mask_file in exclude_list:
                continue
            mask_file = os.path.join(mask_folder, mask_file)
            mask = np.fromfile(mask_file, dtype=bool)
            mask = np.reshape(mask, shape_org)
            OAR_mask = np.logical_or(OAR_mask, mask)
        OAR_mask = np.uint8(OAR_mask * 255)

        output_folder = os.path.join(data_folder, 'patient{}_OAR_paper'.format(i+1))
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        for j in range(shape_org[2]):
            output_file_name = os.path.join(output_folder, '{:03d}.png'.format(j+1))
            plt.imsave(output_file_name, OAR_mask[:, :, j])
        print('patient{} done!'.format(i+1))


def calc_distance():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    num_angles = 20

    angles_mutual = []
    angles_annealing = []
    angles_BOO = []
    for i in range(num_patients):
        patient_E2E = os.path.join(parent_folder, 'patient{}_E2E'.format(i+1))
        beam_angles_annealing_file = os.path.join(patient_E2E, 'beam_angles_annealing_correct.txt')
        beam_angles_BOO_file = os.path.join(patient_E2E, 'beam_angles_E2E.txt')
        if i == 0:
            beam_angles_annealing_file = os.path.join(patient_E2E, 'beam_angles_annealing_correct_additional.txt')
        beam_angles_annealing = read_angles(beam_angles_annealing_file)
        beam_angles_BOO = read_angles(beam_angles_BOO_file)
        vectors_annealing = angle_to_vector(beam_angles_annealing)
        vectors_BOO = angle_to_vector(beam_angles_BOO)

        angular_distance = calc_angular_distance(vectors_annealing, vectors_BOO)
        min_angular_distance0 = np.min(angular_distance, axis=1)
        min_angular_distance1 = np.min(angular_distance, axis=0)
        angle = (np.sum(min_angular_distance0) + np.sum(min_angular_distance1)) / (2 * num_angles)
        angles_mutual.append(angle)

        angular_distance = calc_angular_distance(vectors_annealing, vectors_annealing)
        angle = 0
        for j in range(num_angles - 1):
            for k in range(j+1, num_angles):
                angle += angular_distance[j, k]
        angle /= num_angles * (num_angles - 1) / 2
        angles_annealing.append(angle)

        angular_distance = calc_angular_distance(vectors_BOO, vectors_BOO)
        angle = 0
        for j in range(num_angles - 1):
            for k in range(j + 1, num_angles):
                angle += angular_distance[j, k]
        angle /= num_angles * (num_angles - 1) / 2
        angles_BOO.append(angle)

    for i in range(num_patients):
        print('for patient {}, mutual distance is {} rad, or {} degree'.format(
            i+1, angles_mutual[i], angles_mutual[i] * 180 / np.pi))
    print('\n')
    for i in range(num_patients):
        print('for patient {}, intra distance is {} rad, or {} degree'.format(
            i+1, angles_annealing[i], angles_annealing[i] * 180 / np.pi))
    print('\n')
    for i in range(num_patients):
        print('for patient {}, intra distance is {} rad, or {} degree'.format(
            i + 1, angles_BOO[i], angles_BOO[i] * 180 / np.pi))


def read_angles(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    num_angles = 20
    matrix = np.zeros((num_angles, 2), dtype=np.float32)
    for i in range(num_angles):
        line = lines[i]
        line = line.split(',')
        matrix[i, 0], matrix[i, 1] = float(line[0]), float(line[1])
    return matrix


def angle_to_vector(angles):
    num_angles = 20
    vectors = np.zeros((num_angles, 3), dtype=np.float32)
    vectors[:, 0] = np.sin(angles[:, 0]) * np.cos(angles[:, 1])
    vectors[:, 1] = np.sin(angles[:, 0]) * np.sin(angles[:, 1])
    vectors[:, 2] = np.cos(angles[:, 0])
    return vectors


def calc_angular_distance(vectors0, vectors1):
    vectors0 = np.expand_dims(vectors0, axis=1)
    vectors1 = np.expand_dims(vectors1, axis=0)
    middle = vectors0 * vectors1
    middle = np.sum(middle, axis=2)
    middle = np.arccos(middle)
    return middle


def fuseContourImage():
    num_patients = 6
    image_range = [(50, 1300), (600, 1500)]
    image_size = (image_range[0][1] - image_range[0][0], image_range[1][1] - image_range[1][0])
    # image_size = (750, 450)
    image_size_fuse = (image_size[0]*3, image_size[1]*4, 3)
    fuse = np.zeros(image_size_fuse, dtype=np.float32)

    coords_start_annealing = []
    for i in range(3):
        for j in range(2):
            coord_start = (i * image_size[0], j * image_size[1] * 2)
            coords_start_annealing.append(coord_start)
    coords_start_BOO = [(a[0], a[1]+image_size[1]) for a in coords_start_annealing]

    for i in range(num_patients):
        patientContourAnnealingPath = 'patient{}ContourAnnealing.png'.format(i+1)
        patientContourBOOPath = 'patient{}ContourBOO.png'.format(i+1)
        imageAnnealing_ = plt.imread(patientContourAnnealingPath)
        imageBOO_ = plt.imread(patientContourBOOPath)
        # plt.imshow(imageAnnealing_)
        # plt.show()
        # break
        imageAnnealing = imageAnnealing_[image_range[0][0]: image_range[0][1],
                         image_range[1][0]: image_range[1][1]]
        imageBOO = imageBOO_[image_range[0][0]: image_range[0][1],
                   image_range[1][0]: image_range[1][1]]
        coord_start_annealing = coords_start_annealing[i]
        coord_start_BOO = coords_start_BOO[i]
        fuse[coord_start_annealing[0]:coord_start_annealing[0]+image_size[0],
            coord_start_annealing[1]:coord_start_annealing[1]+image_size[1]] = imageAnnealing
        fuse[coord_start_BOO[0]:coord_start_BOO[0]+image_size[0],
            coord_start_BOO[1]:coord_start_BOO[1]+image_size[1]] = imageBOO

    text_offset_x = 100
    text_offset_y = 50
    text_offset_z = 200
    font = cv2.FONT_HERSHEY_SIMPLEX
    fontSize = 3
    color = (0, 0, 0)
    thickness = 8
    for i in range(num_patients):
        coord_start_annealing = coords_start_annealing[i]
        org = (coord_start_annealing[1] + text_offset_y, coord_start_annealing[0] + text_offset_x)
        text = 'patient {}'.format(i+1)
        fuse = cv2.putText(fuse, text, org, font, fontSize, color, thickness)

        org = (coord_start_annealing[1] + text_offset_y, coord_start_annealing[0] + text_offset_z)
        text = 'annealing'
        fuse = cv2.putText(fuse, text, org, font, fontSize, color, thickness)

    for i in range(num_patients):
        coord_start_BOO = coords_start_BOO[i]
        org = (coord_start_BOO[1] + text_offset_y, coord_start_BOO[0] + text_offset_x)
        text = 'patient {}'.format(i+1)
        fuse = cv2.putText(fuse, text, org, font, fontSize, color, thickness)

        org = (coord_start_BOO[1] + text_offset_y, coord_start_BOO[0] + text_offset_z)
        text = 'baseline'
        fuse = cv2.putText(fuse, text, org, font, fontSize, color, thickness)

    output_file = './contour.png'
    plt.imsave(output_file, fuse)


def draw_dose_axial():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/'
    num_patients = 6
    voxel_size = 2
    append_module = 8
    slice = 'sagittal'  # 'coronal', 'axial' ,
    exclude_list = ['BODY.dat', 'Skin.dat', 'RingStructure.dat']
    color_list = [(255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 255, 0), (212, 255, 127),
                  (87, 207, 227), (0, 255, 127), (30, 105, 210), (237, 149, 100), (34, 139, 34)]
    # colormap = cv2.COLORMAP_SUMMER
    colormap = cv2.COLORMAP_JET
    color_range = 0.8
    mix_weight = 0.8

    output_folder = '/data/qifan/projects_qlyu/EndtoEnd3/papergraph/doseSlices'
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    for i in range(num_patients):
        patient_E2E = os.path.join(parent_folder, 'patient{}_E2E'.format(i+1))
        CT_file = os.path.join(patient_E2E, 'CT.dat')
        shape_file = os.path.join(patient_E2E, 'shape.txt')
        isocenter_file = os.path.join(patient_E2E, 'isocenter.txt')

        # read shape
        with open(shape_file, 'r') as f:
            line = f.readline()
        line = line.split(' ')
        shape = (int(line[0]), int(line[1]), int(line[2]))
        shape_append = (shape[0], shape[1], int(append_module * np.ceil(shape[2] / append_module)))

        # read isocenter
        with open(isocenter_file, 'r') as f:
            line = f.readline()
        line = line.split(',')
        isocenter = (float(line[0]), float(line[1]), float(line[2]))
        # print(shape)
        # print(isocenter)
        # print('\n')
        # calculate central slice idx
        central_slice = (int(isocenter[0]/voxel_size),
                         int(isocenter[1]/voxel_size), int(isocenter[2]/voxel_size))

        # load CT
        CT_ = np.fromfile(CT_file, dtype=np.float32)
        CT = np.zeros(shape_append, dtype=np.float32)
        CT[:, :, :shape[2]] = np.reshape(CT_, shape)
        # load annealing dose and BOO dose
        annealingDoseFile = os.path.join(parent_folder,
            'patient{}_annealing_correct_init'.format(i+1), 'totalDose.dat')
        BOODoseFile = os.path.join(parent_folder,
            'patient{}_BOO_correct'.format(i+1), 'totalDose.dat')
        if i == 0:
            annealingDoseFile = os.path.join(parent_folder,
                'patient1_annealing_correct_additional_init', 'totalDose.dat')
        annealingDose = np.fromfile(annealingDoseFile, dtype=np.float32)
        annealingDose = np.reshape(annealingDose, shape_append)
        BOODose = np.fromfile(BOODoseFile, dtype=np.float32)
        BOODose = np.reshape(BOODose, shape_append)

        # load masks
        mask_folder = os.path.join(parent_folder, 'patient{}_compare'.format(i+1), 'masks')
        files = os.listdir(mask_folder)
        masks = {}
        for file in files:
            if file in exclude_list:
                continue
            fullFile = os.path.join(mask_folder, file)
            mask_ = np.fromfile(fullFile, dtype=bool)
            mask = np.zeros(shape_append, dtype=bool)
            mask[:, :, :shape[2]] = np.reshape(mask_, shape)

            name = file.split('.')[0]
            masks[name] = mask

        # load body mask
        body_mask_file = os.path.join(mask_folder, 'BODY.dat')
        bodyMask_ = np.fromfile(body_mask_file, dtype=bool)
        bodyMask = np.zeros(shape_append, dtype=bool)
        bodyMask[:, :, :shape[2]] = np.reshape(bodyMask_, shape)

        mask_slices = {}
        if slice == 'axial':
            CTSlice = CT[:, :, central_slice[2]]
            annealingDoseSlice = annealingDose[:, :, central_slice[2]]
            BOODoseSlice = BOODose[:, :, central_slice[2]]
            for name, mask in masks.items():
                mask_slices[name] = mask[:, :, central_slice[2]]
            bodyMaskSlice = bodyMask[:, :, central_slice[2]]
        elif slice == 'coronal':
            CTSlice = CT[central_slice[0], :, :]
            annealingDoseSlice = annealingDose[central_slice[0], :, :]
            BOODoseSlice = BOODose[central_slice[0], :, :]
            for name, mask in masks.items():
                mask_slices[name] = mask[central_slice[0], :, :]
            bodyMaskSlice = bodyMask[central_slice[0], :, :]
        elif slice == 'sagittal':
            CTSlice = CT[:, central_slice[1], :]
            annealingDoseSlice = annealingDose[:, central_slice[1], :]
            BOODoseSlice = BOODose[:, central_slice[1], :]
            for name, mask in masks.items():
                mask_slices[name] = mask[:, central_slice[1], :]
            bodyMaskSlice = bodyMask[:, central_slice[1], :]
        else:
            raise 'slice == {}, not correct'.format(slice)

        # calculate contour
        result_ = np.uint8(255 * CTSlice / np.max(CTSlice))
        result_shape = (result_.shape[0], result_.shape[1], 3)
        result = np.zeros(result_shape, dtype=np.uint8)
        result[:, :, 0] = result_
        result[:, :, 1] = result_
        result[:, :, 2] = result_
        color_idx = 0
        for name, mask in mask_slices.items():
            # print(np.max(mask), mask.dtype)
            mask = np.uint8(255*mask)
            ret, thresh = cv2.threshold(mask, 128, 255, 0)
            contours, hierachy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
            cv2.drawContours(result, contours, -1, color_list[color_idx], 1)
            color_idx += 1

        annealingDoseSlice_normalized = np.uint8(255 * annealingDoseSlice / np.max(annealingDoseSlice))
        annealingDoseSlice_normalized = 255 - annealingDoseSlice_normalized
        annealingDoseSlice_normalized = np.uint8(annealingDoseSlice_normalized * color_range + 255*(1-color_range))
        annealingColormap = cv2.applyColorMap(annealingDoseSlice_normalized, colormap)
        annealingColormap[np.logical_not(bodyMaskSlice)] = 0
        BOODoseSlice_normalized = np.uint8(255 * BOODoseSlice / np.max(BOODoseSlice))
        BOODoseSlice_normalized = 255 - BOODoseSlice_normalized
        BOODoseSlice_normalized = np.uint8(BOODoseSlice_normalized * color_range + 255*(1-color_range))
        BOOColormap = cv2.applyColorMap(BOODoseSlice_normalized, colormap)
        BOOColormap[np.logical_not(bodyMaskSlice)] = 0

        annealing_mix = np.uint8(mix_weight * result + (1-mix_weight)*annealingColormap)
        BOO_mix = np.uint8(mix_weight * result + (1-mix_weight)*BOOColormap)

        output_file_annealing = os.path.join(output_folder, 'patient{}Annealing{}.png'.format(i+1, slice))
        output_file_BOO = os.path.join(output_folder, 'patient{}BOO{}.png'.format(i+1, slice))
        plt.imsave(output_file_annealing, annealing_mix)
        plt.imsave(output_file_BOO, BOO_mix)


def combine_dose_image():
    # we only include patient 1 and patient 5
    doseSliceFolder = '/data/qifan/projects_qlyu/EndtoEnd3/papergraph/doseSlices'
    patients = ['patient1', 'patient5']
    methods = ['Annealing', 'BOO']
    slices = ['axial', 'coronal', 'sagittal']
    images = []
    for patient in patients:
        for method in methods:
            for slice in slices:
                file = os.path.join(doseSliceFolder, patient + method + slice + '.png')
                image = plt.imread(file)
                # image = cv2.imread(file)
                images.append(image)
    # for image in images:
    #     print(image.shape, image.dtype)
    #     plt.imshow(image)
    #     plt.show()

    rows = 4
    columns = 3
    size = 200
    fuse_size = (size * rows, size * columns, images[0].shape[2])
    fuse = np.zeros(fuse_size, images[0].dtype)
    canvas = np.zeros((size, size, images[0].shape[2]), dtype=images[0].dtype)
    canvas[:, :, 3] = 1.0
    for row in range(rows):
        for col in range(columns):
            idx = row * columns + col
            image = images[idx]
            if col == 1:
                image = cv2.rotate(image, cv2.ROTATE_90_CLOCKWISE)
            if col == 2:
                image = cv2.rotate(image, cv2.ROTATE_90_CLOCKWISE)
                image = cv2.flip(image, 1)

            margin0 = int((size - image.shape[0]) / 2)
            margin1 = int((size - image.shape[1]) / 2)
            canvas[:, :, :3] = 0
            canvas[margin0: margin0+image.shape[0], margin1: margin1+image.shape[1]] = image
            fuse[row*size: (row+1)*size, col*size: (col+1)*size] = canvas

    offset0 = 10
    offset1 = 20
    offset2 = 40
    for row in range(rows):
        if row in (0, 1):
            patient_name = 'patient 1'
        elif row in (2, 3):
            patient_name = 'patient 5'

        if row in (0, 2):
            method = 'annealing'
        elif row in (1, 3):
            method = 'baseline'

        line = '{}\n{}'.format(patient_name, method)
        plt.text(offset0, offset2 + row * size, line, color='white', size=6)

    output_path = '/data/qifan/projects_qlyu/EndtoEnd3/papergraph/doseSlices.png'
    plt.axis('off')
    plt.imshow(fuse)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)


def draw_water_dose():
    QX_water_out = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/QX_water_out'
    Ryan_water_out = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/Ryan_water_out'
    QX_shape = (200, 200, 200)
    Ryan_shape = (160, 160, 160)
    QX_central = 100
    Ryan_central = 80

    target_size = (200, 200)
    num_beams = 5
    color_range = 0.8
    images = []
    colormap = cv2.COLORMAP_JET

    for i in range(num_beams):
        QX_water_file = os.path.join(QX_water_out, 'waterDose_1.57_{:.2f}.dat'.format(2*np.pi*i/num_beams))
        Ryan_water_file = os.path.join(Ryan_water_out, 'full_matrix_{}.dat'.format((5-i)%5))
        QX_dose = np.fromfile(QX_water_file, dtype=np.float32)
        QX_dose = np.reshape(QX_dose, QX_shape)
        Ryan_dose = np.fromfile(Ryan_water_file, dtype=np.float64)
        Ryan_dose = np.reshape(Ryan_dose, Ryan_shape)
        QX_central_slice = QX_dose[:, :, QX_central].copy()
        Ryan_central_slice = Ryan_dose[Ryan_central, :, :]
        Ryan_central_slice = cv2.resize(Ryan_central_slice, target_size)

        QX_central_slice[:50, :] = 0
        Ryan_central_slice[:50, :] = 0
        QX_central_slice = np.uint8(255 * QX_central_slice / np.max(QX_central_slice))
        Ryan_central_slice = np.uint8(255 * Ryan_central_slice / np.max(Ryan_central_slice))
        QX_slice = np.zeros((target_size[0], target_size[1], 3), dtype=np.uint8)
        QX_slice[:, :, 0] = QX_central_slice
        QX_slice[:, :, 1] = QX_central_slice
        QX_slice[:, :, 2] = QX_central_slice

        Ryan_slice = np.zeros((target_size[0], target_size[1], 3), dtype=np.uint8)
        Ryan_slice[:, :, 0] = Ryan_central_slice
        Ryan_slice[:, :, 1] = Ryan_central_slice
        Ryan_slice[:, :, 2] = Ryan_central_slice

        QX_slice = 255 - QX_slice
        QX_slice = np.uint8(color_range * QX_slice + (1 - color_range) * 255)
        QX_slice = cv2.applyColorMap(QX_slice, colormap)

        Ryan_slice = 255  -Ryan_slice
        Ryan_slice = np.uint8(color_range * Ryan_slice + (1 - color_range) * 255)
        Ryan_slice = cv2.applyColorMap(Ryan_slice, colormap)

        images.append(QX_slice)
        images.append(Ryan_slice)

    rows = 3
    columns = 4
    size = 200
    fuse_size = (rows*size, columns*size, 3)
    fuse = np.ones(fuse_size, dtype=np.uint8)
    fuse *= 255
    for row in range(rows):
        for col in range(columns):
            idx = row * columns + col
            if idx == 10:
                break
            image = images[idx]
            fuse[row*size:(row+1)*size, col*size:(col+1)*size, :] = image

    # text
    offset0 = 30
    offset1 = 10
    for row in range(rows):
        for col in range(columns):
            idx = row * columns + col
            if idx % 2 == 0:
                text = 'FCBB'
            elif idx % 2 == 1:
                text = 'CCCS'
            plt.text(offset1 + col * size, offset0 + row * size, text, color='white', size=8)
    output_file = '/data/qifan/projects_qlyu/EndtoEnd3/papergraph/waterDose.png'
    plt.axis('off')
    plt.imshow(fuse)
    plt.savefig(output_file, bbox_inches='tight', dpi=300)


def plot_water_dose():
    QX_water_file = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/QX_water_out/waterDose_1.57_0.00.dat'
    Ryan_water_file = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/Ryan_water_out/full_matrix_0.dat'
    QX_shape = (200, 200, 200)
    Ryan_shape = (160, 160, 160)

    QX_dose = np.fromfile(QX_water_file, dtype=np.float32)
    QX_dose = np.reshape(QX_dose, QX_shape)
    Ryan_dose = np.fromfile(Ryan_water_file, dtype=np.float64)
    Ryan_dose = np.reshape(Ryan_dose, Ryan_shape)
    QX_slice = QX_dose[:, :, 100]
    Ryan_slice = Ryan_dose[80, :, :]

    size = 200
    shape = (200, 200)
    Ryan_slice = cv2.resize(Ryan_slice, shape)
    QX_slice[:50, :] = 0
    Ryan_slice[:50, :] = 0

    # align
    A = np.sum(Ryan_slice * Ryan_slice)
    B = np.sum(Ryan_slice * QX_slice)
    k = B / A
    Ryan_slice *= k

    fig, axes = plt.subplots(1, 2)

    QX_axial = QX_slice[50:, 100].copy()
    Ryan_axial = Ryan_slice[50:, 100].copy()
    points = 150
    roof = np.max(QX_axial)
    QX_axial /= roof
    Ryan_axial /= roof
    axes[0].plot(np.arange(points)*0.2, QX_axial, color='tab:blue', linestyle='-')
    axes[0].plot(np.arange(points)*0.2, Ryan_axial, color='tab:blue', linestyle='--')
    axes[0].set_xlabel('depth/cm', fontsize=16)
    axes[0].set_ylabel('dose (a.u.)', fontsize=16)
    axes[0].legend(['FCBB', 'CCCS'], fontsize=14)

    axes[0].set_title('axial dose profile'.format(1), fontsize=20)
    axes[0].tick_params(axis='both', which='major', labelsize=14)
    axes[0].tick_params(axis='both', which='major', labelsize=14)

    depths = [2, 5, 10, 15, 20]
    labels = ['{}cm'.format(a) for a in depths]
    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    depths = [50 + a * 5 for a in depths]

    # QX transverse
    FOV = 80
    margin = int((size - FOV) / 2)
    for i, depth in enumerate(depths):
        line = QX_slice[depth, :]
        line = line[margin: margin+FOV]
        axes[1].plot((np.arange(FOV)-FOV/2)*0.2, line, color=color_list[i], linestyle='-')
    axes[1].legend(labels, fontsize=14)

    for i, depth in enumerate(depths):
        line = Ryan_slice[depth, :]
        line = line[margin: margin+FOV]
        axes[1].plot((np.arange(FOV) - FOV / 2 - 0.5) * 0.2, line, color=color_list[i], linestyle='--')
        axes[1].plot()
    axes[1].set_title('transverse dose profile', fontsize=20)
    axes[1].set_xlabel('offset/cm', fontsize=16)
    axes[1].tick_params(axis='both', which='major', labelsize=14)
    axes[1].tick_params(axis='both', which='major', labelsize=14)

    fig.set_figheight(5)
    fig.set_figwidth(13)
    fig.tight_layout()
    outputFile = '/data/qifan/projects_qlyu/EndtoEnd3/papergraph/waterDoesPlot.pdf'
    fig.savefig(outputFile)

    # plt.imshow(QX_slice)
    # plt.show()
    # plt.imshow(Ryan_slice)
    # plt.show()


def draw_patient_dose():
    Ryan_template = '/data/qifan/projects_qlyu/dose-calculation/patients/patient{}/dosecalc'
    QX_template = '/data/qifan/projects_qlyu/dose-calculation/patients/patient{}_E2E'
    shapes_org = [(200, 200, 197), (200, 200, 138), (171, 171, 112),
                  (170, 170, 113), (200, 200, 152), (260, 260, 128)]
    append_module = 8
    shapes_append = [(a[0], a[1], append_module*int(np.ceil(a[2]/append_module))) for a in shapes_org]
    shapes_Ryan = [(160, 160, 158), (160, 160, 111), (137, 137, 90),
                   (136, 136, 91), (160, 160, 122), (208, 208, 103)]
    num_patients = 6
    beams_per_patient = 5
    for i in range(num_patients):
        Ryan_folder = Ryan_template.format(i+1)
        QX_folder = QX_template.format(i+1)
        shape_Ryan = shapes_Ryan[i]
        shape_Ryan = (shape_Ryan[2], shape_Ryan[0], shape_Ryan[1])
        shape_QX = shapes_append[i]
        shape_org = shapes_org[i]

        central_slice_idx_QX = int(shape_org[2]/2)
        central_slice_idx_Ryan = int(shape_Ryan[0] / 2)
        for j in range(1, beams_per_patient):
            Ryan_beam = os.path.join(Ryan_folder, 'beam{}.dat'.format((beams_per_patient-j)%beams_per_patient+1))
            QX_beam = os.path.join(QX_folder, 'waterDose_1.57_{:.2f}.dat'.format(np.pi*2*j/beams_per_patient))
            Ryan_dose = np.fromfile(Ryan_beam, dtype=np.float64)
            Ryan_dose = np.reshape(Ryan_dose, shape_Ryan)
            Ryan_dose = np.transpose(Ryan_dose, (1, 2, 0))
            Ryan_dose = np.flip(Ryan_dose, axis=2)
            Ryan_dose = trilinear_interpolate(Ryan_dose, shape_org)

            output_file = os.path.join(Ryan_folder, 'beam{}_resize.dat'.format(
                (beams_per_patient-j)%beams_per_patient+1))
            Ryan_dose.tofile(output_file)
            print(i, j)


def trilinear_interpolate(input, new_shape):
    output_array = np.zeros(new_shape, dtype=input.dtype)
    old_shape = input.shape
    for i in range(new_shape[0]):
        I = (i + 0.5) / new_shape[0] * old_shape[0] - 0.5
        for j in range(new_shape[1]):
            J = (j + 0.5) / new_shape[1] * old_shape[1] - 0.5
            for k in range(new_shape[2]):
                K = (k + 0.5) / new_shape[2] * old_shape[2] - 0.5
                value = 0
                for l in range(2):
                    IDX = int(np.floor(I)) + l
                    coeffi = 1 - abs(I - IDX)
                    if IDX < 0 or IDX >= old_shape[0]:
                        continue
                    for m in range(2):
                        JDX = int(np.floor(J)) + m
                        coeffj = 1 - abs(J - JDX)
                        if JDX < 0 or JDX >= old_shape[1]:
                            continue
                        for n in range(2):
                            KDX = int(np.floor(K)) + n
                            coeffk = 1 - abs(K - KDX)
                            if KDX < 0 or KDX >= old_shape[2]:
                                continue
                            value += coeffi * coeffj * coeffk * input[IDX, JDX, KDX]
                output_array[i, j, k] = value
    return output_array


def draw_patient_dose1():
    Ryan_template = '/data/qifan/projects_qlyu/dose-calculation/patients/patient{}/dosecalc'
    QX_template = '/data/qifan/projects_qlyu/dose-calculation/patients/patient{}_E2E'
    body_mask_template = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient{}_compare/masks/BODY.dat'
    CT_template = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient{}_E2E/CT.dat'
    shapes_org = [(200, 200, 197), (200, 200, 138), (171, 171, 112),
                  (170, 170, 113), (200, 200, 152), (260, 260, 128)]
    append_module = 8
    shapes_append = [(a[0], a[1], append_module * int(np.ceil(a[2] / append_module))) for a in shapes_org]
    num_patients = 6
    beams_per_patient = 5
    # for i in range(num_patients):
    images_QX = []
    images_Ryan = []
    colormap = cv2.COLORMAP_JET
    color_range = 0.8
    mix_weight = 0.8
    for i in (0, 4):
        Ryan_folder = Ryan_template.format(i+1)
        QX_folder = QX_template.format(i+1)
        shape_org = shapes_org[i]
        shape_append = shapes_append[i]
        central_slice = int(shape_org[2] / 2)

        mask_file = body_mask_template.format(i+1)
        mask = np.fromfile(mask_file, dtype=bool)
        mask = np.reshape(mask, shape_org)
        mask_slice = mask[:, :, central_slice]
        mask_slice = np.expand_dims(mask_slice, axis=2)

        CT_file = CT_template.format(i+1)
        CT = np.fromfile(CT_file, dtype=np.float32)
        CT = np.reshape(CT, shape_org)
        CT_slice_ = CT[:, :, central_slice].copy()
        CT_slice_ = CT_slice_ / np.max(CT_slice_)
        CT_slice_ = np.uint8(255 * CT_slice_)
        CT_slice = np.zeros((CT_slice_.shape[0], CT_slice_.shape[1], 3), dtype=np.uint8)
        CT_slice[:, :, 0] = CT_slice_
        CT_slice[:, :, 1] = CT_slice_
        CT_slice[:, :, 2] = CT_slice_

        for j in range(beams_per_patient):
            Ryan_file = os.path.join(Ryan_folder, 'beam{}_resize.dat'.format((beams_per_patient-j) % 5+1))
            Ryan_dose = np.fromfile(Ryan_file, dtype=np.float64)
            Ryan_dose = np.reshape(Ryan_dose, shape_org)
            Ryan_dose *= mask

            QX_file = os.path.join(QX_folder, 'waterDose_1.57_{:.2f}.dat'.format(np.pi*2*j/beams_per_patient))
            QX_dose = np.fromfile(QX_file, dtype=np.float32)
            QX_dose = np.reshape(QX_dose, shape_append)
            QX_dose = QX_dose[:, :, :shape_org[2]]
            QX_dose *= mask
            roof = np.max(QX_dose)
            QX_dose /= roof

            # normalize
            k = np.sum(Ryan_dose * QX_dose) / np.sum(Ryan_dose * Ryan_dose)
            Ryan_dose *= k

            central_slice_Ryan = Ryan_dose[:, :, central_slice].copy()
            central_slice_Ryan = color_range * (1 - central_slice_Ryan) + (1 - color_range)
            # central_slice_Ryan[central_slice_Ryan>1] = 1
            # central_slice_Ryan[central_slice_Ryan<0] = 0
            central_slice_Ryan = np.uint8(central_slice_Ryan * 255)
            # print(np.max(central_slice_Ryan), np.min(central_slice_Ryan))
            colormap_Ryan = cv2.applyColorMap(central_slice_Ryan, colormap)
            colormap_Ryan = colormap_Ryan * mask_slice
            colormap_Ryan = mix_weight * CT_slice + (1 - mix_weight) * colormap_Ryan
            colormap_Ryan = np.uint8(colormap_Ryan)

            central_slice_QX = QX_dose[:, :, central_slice].copy()
            central_slice_QX = color_range * (1 - central_slice_QX) + (1 - color_range)
            central_slice_QX = np.uint8(central_slice_QX * 255)
            # print(np.max(central_slice_QX), np.min(central_slice_QX))
            colormap_QX = cv2.applyColorMap(central_slice_QX, colormap)
            colormap_QX = colormap_QX * mask_slice
            colormap_QX = mix_weight * CT_slice + (1 - mix_weight) * colormap_QX
            colormap_QX = np.uint8(colormap_QX)

            images_QX.append(colormap_QX)
            images_Ryan.append(colormap_Ryan)

    size = 200
    cols = 4
    rows = 5
    canvas_size = (rows*size, cols*size, 3)
    canvas = np.zeros(canvas_size, dtype=np.uint8)

    # QX
    columns = 2
    for row in range(rows):
        for col in range(columns):
            idx = col * rows + row
            column = col * 2
            # print(images_QX[idx].shape)
            canvas[row*size:(row+1)*size, column*size:(column+1)*size] = images_QX[idx]

    # Ryan
    for row in range(rows):
        for col in range(columns):
            idx = col * rows + row
            column = col * 2 + 1
            # print(images_Ryan[idx].shape)
            canvas[row*size:(row+1)*size, column*size: (column+1)*size] = images_Ryan[idx]

    legends = ['patient 1 FCBB', 'patient 1 CCCS', 'patient 5 FCBB', 'patient 5 CCCS']
    row_offset = 10
    col_offset = 30
    for i in range(4):
        plt.text(row_offset + i * size, col_offset, legends[i], color='white', fontsize=6)
    plt.tight_layout()
    plt.axis('off')
    plt.imshow(canvas)
    output_path = '/data/qifan/projects_qlyu/EndtoEnd3/papergraph/patientDose.png'
    plt.savefig(output_path, bbox_inches='tight', dpi=300)


def calc_dose_diff():
    Ryan_template = '/data/qifan/projects_qlyu/dose-calculation/patients/patient{}/dosecalc'
    QX_template = '/data/qifan/projects_qlyu/dose-calculation/patients/patient{}_E2E'
    body_mask_template = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient{}_compare/masks/BODY.dat'
    CT_template = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient{}_E2E/CT.dat'
    shapes_org = [(200, 200, 197), (200, 200, 138), (171, 171, 112),
                  (170, 170, 113), (200, 200, 152), (260, 260, 128)]
    append_module = 8
    shapes_append = [(a[0], a[1], append_module * int(np.ceil(a[2] / append_module))) for a in shapes_org]
    num_patients = 6
    beams_per_patient = 5
    PSNR_total = 0
    for i in range(num_patients):
        QX_folder = QX_template.format(i+1)
        Ryan_folder = Ryan_template.format(i+1)
        shape_org = shapes_org[i]
        shape_append = shapes_append[i]
        central_slice_0 = int(shape_org[0]/2)
        central_slice_1 = int(shape_org[1]/2)

        mask_file = body_mask_template.format(i+1)
        mask = np.fromfile(mask_file, dtype=bool)
        mask = np.reshape(mask, shape_org)
        for j in range(beams_per_patient):
            QX_file = os.path.join(QX_folder, 'waterDose_1.57_{:.2f}.dat'.format(np.pi * 2 * j/beams_per_patient))
            QX_dose = np.fromfile(QX_file, dtype=np.float32)
            QX_dose = np.reshape(QX_dose, shape_append)
            QX_dose = QX_dose[:, :, :shape_org[2]]
            QX_dose *= mask

            Ryan_file = os.path.join(Ryan_folder, 'beam{}_resize.dat'.format(
                (beams_per_patient - j) % beams_per_patient + 1))
            Ryan_dose = np.fromfile(Ryan_file, dtype=np.float64)
            Ryan_dose = np.reshape(Ryan_dose, shape_org)
            Ryan_dose *= mask

            # # alignment check
            # central_slice_0_QX = QX_dose[central_slice_0, :, :]
            # central_slice_0_Ryan = Ryan_dose[central_slice_0, :, :]
            # central_slice_1_QX = QX_dose[:, central_slice_1, :]
            # central_slice_1_Ryan = Ryan_dose[:, central_slice_1, :]
            # plt.imshow(central_slice_0_QX)
            # plt.show()
            # plt.imshow(central_slice_0_Ryan)
            # plt.show()
            # plt.imshow(central_slice_1_QX)
            # plt.show()
            # plt.imshow(central_slice_1_Ryan)
            # plt.show()

            # match
            k = np.sum(Ryan_dose * QX_dose) / np.sum(Ryan_dose * Ryan_dose)
            Ryan_dose *= k
            diff = QX_dose - Ryan_dose
            squared = np.sum(diff * diff)
            mse = squared / np.sum(mask)
            roof = np.max(QX_dose)
            PSNR = 10 * np.log(roof**2 / mse)
            PSNR_total += PSNR
            print('patient {} beam {}: PSNR {}'.format(i, j, PSNR))
    PSNR_avg = PSNR_total / (num_patients * beams_per_patient)
    print('average PSNR: {}'.format(PSNR_avg))


if __name__ == '__main__':
    # draw_annealing_loss()
    # draw_compare_E2E()
    # draw_compare_BOO()
    # draw_DVH()
    # draw_masks()
    # calc_distance()
    # fuseContourImage()
    # draw_dose_axial()
    # draw_dose_axial()
    # combine_dose_image()
    # draw_water_dose()
    # plot_water_dose()
    # draw_patient_dose()
    # draw_patient_dose1()
    calc_dose_diff()