import sys
import numpy as np
import matplotlib.pyplot as plt
from fipy import *
from scipy.stats import norm
import time


def run_solution(length, dx):
    nx = length // dx
    D = 50.
    steps = 1000
    mesh = Grid1D(nx=nx, dx=dx)
    phi = CellVariable(name="numeric solution", mesh=mesh, value=0.)
    solution = CellVariable(name="analytic solution", mesh=mesh, value=0.)
    x = mesh.cellCenters[0]

    init_profile = 0.3
    phi.setValue(init_profile)
    solution.setValue(init_profile)

    phi.faceGrad.constrain(0, mesh.facesRight)
    phi.constrain(0, mesh.facesLeft)
    
    SRC = 0  # reaction part
    eqX = TransientTerm() == ExplicitDiffusionTerm(coeff=D) + SRC
    eqI = TransientTerm() == ImplicitDiffusionTerm(coeff=D) + SRC
    eqCN = eqX + eqI
    
    # solution timesteps
    ref_ts = 0.9 * dx**2 / (2 * D)
    ts = ref_ts * 5
    fake_time = 0

    real_time = 0
    start, end = 0, 0

    for step in range(steps):
        if step == 0:
            fake_time += ts
            eqCN.solve(var=phi, dt=ts)
            solution.setValue(init_profile * (2 * norm.cdf(x / np.sqrt(2 * D * fake_time)) - 1))
        else:
            fake_time += ts
            start = time.time()
            eqCN.solve(var=phi, dt=ts)
            end = time.time()
            solution.setValue(init_profile * (2 * norm.cdf(x / np.sqrt(2 * D * fake_time)) - 1))
            real_time += end - start
    return real_time


steps_number = np.linspace(1e2, 3e4, 100)
times = np.array([])
with open("time.txt", 'w') as f:
    f.truncate()
    f.write("steps_number time\n")
    for k in steps_number:
        print k
        f.write(str(k) + ' ' + str(run_solution(50, 50 / k)) + '\n')
        # times = np.append(times, run_solutions(50, 50 / k))

# file = open("time.txt", "w")
# file.truncate()
# file.write("first line is number of grid steps, the second is times\n")
# file.write(str(list(steps_number)) + '\n')
# file.write(str(list(times)))
# file.close()

# plt.plot(steps_number[1:], times[1:])
# plt.ylabel("time(sec)")
# plt.xlabel("grid steps number")
# plt.grid()
# plt.show()
