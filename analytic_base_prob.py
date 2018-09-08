import sys
import numpy as np
import matplotlib.pyplot as plt
from fipy import *
from scipy.stats import norm

# print 'gauss:', norm.pdf(0)
# print 'cdf:', norm.cdf(0)

nx = 100
dx = 0.1

mesh = Grid1D(nx=nx, dx=dx)
x = mesh.cellCenters[0]
print norm.cdf(x)

c0 = np.exp(-x) * x
solution = CellVariable(name="main var", mesh=mesh, value=0.)
solution.setValue(c0)

ax = plt.axes()
viewer = Viewer(vars=(solution), datamin=0., datamax=1., axes=ax)
ax.grid()

steps = 100000
ref_ts = 0.9 * dx**2 / 2
ts = ref_ts * 15
time = 0

for step in range(steps):
    time += ts
    solution.setValue(c0 * norm.cdf(x / np.sqrt(time) / 2))
    if step % 20 == 0:
        # print solution.value[499]
        a = c0 * norm.cdf(x / np.sqrt(time) / 2)
        # print a[499]
        viewer.plot()
        print step


