#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
import numpy as np

'''for i in range(0,3):
	x,y,z,dens,vx,vy,vz,num,cpu = np.loadtxt('particles' + str(i), skiprows=5,unpack=True)

	fig = plt.figure()
	ax  = fig.add_subplot(111, projection='3d')

	ax.scatter(x,y,z)
	plt.show()
	plt.close()'''

for i in range(0,10):
        x,y,z,dens,vx,vy,vz,num,cpu = np.loadtxt(('particles0000' + str(i)), skiprows=5,unpack=True)

        fig = plt.figure()
        ax  = fig.add_subplot(111)

        ax.scatter(x,y)
        plt.show()
        plt.close() 
for i in range(10,20):
        x,y,z,dens,vx,vy,vz,num,cpu = np.loadtxt(('particles000' + str(i)), skiprows=5,unpack=True)

        fig = plt.figure()
        ax  = fig.add_subplot(111)

        ax.scatter(x,y)
        plt.show()
        plt.close()
