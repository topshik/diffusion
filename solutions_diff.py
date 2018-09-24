import sys
import numpy as np
import matplotlib.pyplot as plt
from fipy import *
from scipy.stats import norm


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
    
    # timesteps
    ref_ts = 0.9 * dx**2 / (2 * D)
    ts = ref_ts * 5
    time = 0

    plt.figure(figsize=(15, 6))
    ax1 = plt.subplot(2, 2, 1)
    ax1.set_title("numeric solution")
    ax2 = plt.subplot(2, 2, 2)
    ax2.set_title("analytic solution")
    ax3 = plt.subplot(2, 2, 3)
    ax3.set_title("total SE plot")
    ax4 = plt.subplot(2, 2, 4)
    ax4.set_title("SE from coordinate")

    viewer_num = Viewer(vars=(phi), datamin=0., datamax=.5, axes=ax1)
    viewer_ana = Viewer(vars=(solution), datamin=0., datamax=.5, axes=ax2)

    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()

    se_coord_arr = np.array([0])
    total_se_arr = np.array([0])
    time_arr = np.array([0])

    for step in range(steps):
        time += ts
        time_arr = np.append(time_arr, time)
        eqCN.solve(var=phi, dt=ts)
        solution.setValue(init_profile * (2 * norm.cdf(x / np.sqrt(2 * D * time)) - 1))
        se_coord_arr = (np.array(phi) - np.array(solution)) ** 2
        total_se_arr = np.append(total_se_arr, se_coord_arr.sum())
        if step % 20 == 0:
            viewer_num.plot()
            viewer_ana.plot()
            ax3.clear()
            ax3.grid()
            ax3.set_title("total SE plot")
            ax3.plot(time_arr, total_se_arr)

            ax4.cla()
            ax4.grid()
            ax4.set_title("SE from coordinate")
            ax4.plot(np.array(x), se_coord_arr)

            print step


run_solution(100, 0.1)
