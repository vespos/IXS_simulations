import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')



dat = np.loadtxt('./dat/optim.dat', delimiter=',')
xlab = 'pitch correction'
ylab = 'intensity'

fig, ax = plt.subplots()
ax.plot(dat[:,0], dat[:,1])
ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
plt.show()

plt.savefig('./optim.png')