import sys
import numpy as np
import matplotlib.pyplot as plt
from fipy import *
from scipy.stats import norm

nx = 1000
dx = 0.1

# create simple grid at the cut
mesh = Grid1D(nx=nx, dx=dx)
# determine variable at the grid and point initial value
phi = CellVariable(name="numeric var", mesh=mesh, value=0.)
solution = CellVariable(name="analytic var", mesh=mesh, value=0.)

# determine values in cells in respect to cells centers
x = mesh.cellCenters[0]

# set initial profile
# c0 = np.exp(-x) * x
c0 = 0.3
phi.setValue(c0)
solution.setValue(c0)
# you can point values explicitly with phi.value property

D = 50.
steps = 10000

# zero flux. Uncomment it if you want.
phi.faceGrad.constrain(0, mesh.facesRight)

# zero at bounary
phi.constrain(0, mesh.facesLeft)
# phi.constrain(0, mesh.facesRight)


# Defining reaction with diffusion
SRC = 0

# equation with explicit scheme
eqX = TransientTerm() == ExplicitDiffusionTerm(coeff=D) + SRC
# equation with implicit scheme
eqI = TransientTerm() == ImplicitDiffusionTerm(coeff=D) + SRC

# equation with crank-nicolson scheme
eqCN = eqX + eqI

# defining timestep via empiric equation
# you can use less timestep but do not use bigger one
ref_ts = 0.9 * dx**2 / (2 * D)
# if you use CN scheme, use ref_ts*15 as timestep
# if you use explicit scheme use ref_ts as timestep
ts = ref_ts * 5
time = 0

plt.figure(figsize=(15, 6))
ax1 = plt.subplot(1, 2, 1)
viewer_num = Viewer(vars=(phi), datamin=0., datamax=.5, axes=ax1)
ax2 = plt.subplot(1, 2, 2)
viewer_ana = Viewer(vars=(solution), datamin=0., datamax=.5, axes=ax2)
ax1.grid()
ax2.grid()


for step in range(steps):
    time += ts
    eqCN.solve(var=phi, dt=ts)
    solution.setValue(c0 * (2 * norm.cdf(x / np.sqrt(2 * D * time)) - 1))
    if step % 100 == 0:
        viewer_num.plot()
        viewer_ana.plot()
        print step
#    print phi.value
'''

a = []
for step in range(steps):
    time += ts
    a.append(np.array(c0 * (2 * norm.cdf(x / np.sqrt(2 * D * time)) - 1))[50])
print a
plt.ylim(ymin=0, ymax=0.1)
plt.plot(np.array(a))
plt.show()
'''
