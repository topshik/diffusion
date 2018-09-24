import matplotlib.pyplot as plt
import numpy as np

with open("time.txt", 'r') as f:
    line = f.readline()
    data = f.read().split('\n')[:-1]
    steps_number = [float(data[i].split(' ')[0]) for i in range(len(data))]
    times = [float(data[i].split(' ')[1]) for i in range(len(data))]

    print(steps_number)
    plt.plot(steps_number, times, '.')
    plt.ylabel("time, seconds")
    plt.xlabel("number of grid cells")
    plt.grid()
    plt.show()
