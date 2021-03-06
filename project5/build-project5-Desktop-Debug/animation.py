from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

l0 = 1.0/((len(np.loadtxt('data_time.dat')[0])-1)) #---------------!!!!!!!!!!!! if error check this!!
print "l0 =", 1.0/((len(np.loadtxt('data_time.dat')[0])-1))
dt = (l0*l0)/(2.0)
T = 0.5
timesteps = int(T/dt) # update manually!
size = int((1/l0)) + 1


def generate(X, Y, i):

    return np.loadtxt('data_time.dat')[i*(size) -size:i*(size)]

fig = plt.figure()
ax = axes3d.Axes3D(fig)

xs = np.linspace(0, 1, (1/l0)+1)
ys = np.linspace(0, 1, (1/l0)+1)
X, Y = np.meshgrid(xs, ys)
Z = generate(X, Y, 1)
wframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
ax.set_zlim(0,1)



def update(i, ax, fig):
    ax.cla()

    Z = generate(X, Y, i+1) #put in data[i]
    wframe = ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)
    ax.set_zlim(0,1)
    return wframe,

ani = animation.FuncAnimation(fig, update, 
        frames=xrange(timesteps-1), 
        fargs=(ax, fig), interval=100)
plt.show()



#ani.save('explicit_solver.mp4', fps=10, writer="mencoder")



"""

xs = np.linspace(0, 1, (1/l0)+1)
ys = np.linspace(0, 1, (1/l0)+1)

print "xs ", xs
print " data", np.loadtxt('data_time.dat')
print "size ", size

for i in range(0, timesteps):
	print i
	#for j in range(size):
	#	for k in range(size):
	data= np.loadtxt('data_time.dat')[i*(size) -size:i*(size)]
	print "data ", data
			#print "indexes ", i*j, i*k

"""






