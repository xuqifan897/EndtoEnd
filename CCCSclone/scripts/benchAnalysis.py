import numpy as np

def main():
    file = '/data/qifan/projects/EndtoEnd/results/CCCSBench/benchlog.txt'
    with open(file, 'r') as f:
        lines = f.readlines()
    lines = lines[46:]
    total = 0.
    nbeams = 100
    time_array = np.zeros(nbeams)
    for i in range(nbeams):
        line_idx = 2 * i + 1
        line = lines[line_idx]
        segments = line.split(' ')
        while('' in segments):
            segments.remove('')
        time = segments[-2]
        time = float(time)
        time_array[i] = time
    avg = np.average(time_array)
    std = np.std(time_array)
    print("the average dose calculation time for a beam: {} ms, the standard deviation is {} ms".format(avg, std))

if __name__ == '__main__':
    main()