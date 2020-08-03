import numpy as np
import matplotlib.pyplot as plt

#Input parameters
mirrorAngle = 0.0035 # rad
mirrorROC = 500 # m
F = mirrorROC * np.sin(mirrorAngle)/2.0  #focal lenght for grazing incidence mirrors
D = 1.0 # m Distance between mirrors

Z = np.arange(0., 150.0, .1) # distance for plot
pos = np.zeros(len(Z)) # outputs position
ang = np.zeros(len(Z)) # outputs angle

test = np.transpose([0.0,0.001]) #test array 
M0 = [[1,0],[-1.0/F, 1]] #first mirror
M1 = [[1,0],[-1.0/F, 1]] #second mirror
mirrorSeperation = [[1,D],[0,1]]

for i in range(0,len(Z)):
	Z0 = [[1,Z[i]],[0,1]]
	Z1 = [[1,Z[i]],[0,1]]
	Mat = np.dot(Z1,np.dot(M1,np.dot(mirrorSeperation,np.dot(M0,Z0))))
	pos[i],ang[i] = np.dot(Mat,test)


fig, ax = plt.subplots()
ax.plot(Z,pos) #x is distance, y is position 
ax.set_xlabel('Distance [m]')
ax.set_ylabel('drift from axis [m]')
# fig.savefig('position.png')
fig.show()

fig2, ax2 = plt.subplots()
ax2.plot(Z,ang) #x is distance, y is position 
ax2.set_xlabel('Distance [m]')
ax2.set_ylabel('angle [rad]')
# fig2.savefig('angle.png')
fig2.show()