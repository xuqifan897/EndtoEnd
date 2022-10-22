import os
import numpy as np
import matplotlib.pyplot as plt
import mat73


def draw_annealing_loss():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/'
    num_patients = 6
    iterations = 1000
    num_beams = 20
    num_perturbations = 2
    startPoint = 200
    for i in range(num_patients):
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

        plt.plot(np.arange(startPoint, iterations), doseLoss[startPoint: iterations])
        plt.xlabel('iterations')
        plt.ylabel('dose loss (a.u.)')
        plt.ticklabel_format(axis='y', scilimits=[0, 0])
        plt.title('patient {}'.format(i+1))
        plt.show()

    # additional experiment on patient 1
    optimize_annealing_folder = os.path.join(global_folder, 'patient11_annealing_correct')
    doseLossFile = os.path.join(optimize_annealing_folder, 'DoseLoss.dat')
    doseLoss = np.fromfile(doseLossFile, dtype=np.float32)
    plt.plot(np.arange(startPoint, iterations), doseLoss[startPoint: iterations])
    plt.xlabel('iterations')
    plt.ylabel('dose loss (a.u.)')
    plt.ticklabel_format(axis='y', scilimits=[0, 0])
    plt.title('patient {}'.format(1))
    plt.show()


def draw_compare_E2E():
    global_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    iterations = 200
    num_beams = 20
    startPoint = 50
    flag = 'dose'
    # flag = 'total'
    roof_scale = 4

    for i in range(num_patients):
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
            plt.plot(np.arange(iterations), BOODoseLoss)
            plt.plot(np.arange(iterations), annealingDoseLoss)
            plt.ylim([0, roof])
            plt.legend(['BOO', 'annealing'])
            plt.xlabel('iterations')
            plt.ylabel('dose loss (a.u.)')
            plt.title('patient {}'.format(i+1))
            plt.show()
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
        plt.plot(np.arange(iterations), BOODoseLoss)
        plt.plot(np.arange(iterations), annealingDoseLoss)
        plt.ylim([0, roof])
        plt.legend(['BOO', 'annealing'])
        plt.xlabel('iterations')
        plt.ylabel('dose loss (a.u.)')
        plt.title('patient 1')
        plt.show()
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


def draw_compare_BOO():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation'
    num_patients = 6
    iterations = 2000
    start = 20
    stop = 100
    flag = 'dose'
    # flag = 'total'

    for i in range(num_patients):
        patient_folder = os.path.join(parent_folder, 'patient{}_compare'.format(i+1))
        BOO_file = os.path.join(patient_folder, 'optimize_BOO', 'polish_result.mat')
        annealing_file = os.path.join(patient_folder, 'optimize_annealing', 'polish_result.mat')
        BOO_mat = mat73.loadmat(BOO_file)
        annealing_mat = mat73.loadmat(annealing_file)

        # plot costsDF, i.e., dose loss
        if flag == 'dose':
            BOO_costsDF = BOO_mat['result']['costsDF_polish']
            annealing_costsDF = annealing_mat['result']['costsDF_polish']
            plt.plot(np.arange(start, stop), BOO_costsDF[start:stop])
            plt.plot(np.arange(start, stop), annealing_costsDF[start:stop])
            plt.legend(['BOO', 'annealing'])
            plt.xlabel('iterations')
            plt.ylabel('dose loss (a.u.)')
            plt.ticklabel_format(axis='y', scilimits=[0, 0])
            plt.title('patient {}'.format(i+1))
            plt.show()

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


if __name__ == '__main__':
    # draw_annealing_loss()
    # draw_compare_E2E()
    draw_compare_BOO()