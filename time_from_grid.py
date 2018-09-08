import sys
import numpy as np
import matplotlib.pyplot as plt
from fipy import *
from scipy.stats import norm
import time


def run_solutions(length, dx):
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

    # plt.figure(figsize=(15, 6))
    # ax1 = plt.subplot(1, 2, 1)
    # ax1.set_title("numeric solution")
    # ax2 = plt.subplot(1, 2, 2)
    # ax2.set_title("analytic solution")

    # viewer_num = Viewer(vars=(phi), datamin=0., datamax=.5, axes=ax1)
    # viewer_ana = Viewer(vars=(solution), datamin=0., datamax=.5, axes=ax2)

    # ax1.grid()
    # ax2.grid()

    real_time = 0
    start, end = 0, 0

    for step in range(steps):
        if step == 0:
            fake_time += ts
            eqCN.solve(var=phi, dt=ts)
            solution.setValue(init_profile * (2 * norm.cdf(x / np.sqrt(2 * D * fake_time)) - 1))
            # if step % 20 == 0:
                # viewer_num.plot()
                # viewer_ana.plot()
                # print step
        else:
            fake_time += ts
            start = time.time()
            eqCN.solve(var=phi, dt=ts)
            solution.setValue(init_profile * (2 * norm.cdf(x / np.sqrt(2 * D * fake_time)) - 1))
            end = time.time()
            real_time += end - start
            # if step % 100 == 0:
                # viewer_num.plot()
                # viewer_ana.plot()
                # print step
    return real_time

grid_steps = np.linspace(0.005, 0.5, 100)
times = np.array([])
i = 0
for grid_step in grid_steps:
    print i
    times = np.append(times, run_solutions(50, grid_step))
    i += 1

file = open("time.txt", "w")
file.write("frist line is grid steps, the second is times\n")
file.write(str(list(grid_steps)) + '\n')
file.write(str(list(times)))
file.close()

plt.plot(grid_steps[1:], times[1:])
plt.ylabel("time(sec)")
plt.xlabel("grid step size")
plt.grid()
plt.show()
