import matplotlib.pyplot as plt
from math import *
from numpy import*

def readfile(name):
	data = []
	infile = open('%s.dat' % name, 'r')
	numberoflines = MagicNumbers(name)[0]
	numberofobjects = MagicNumbers(name)[1]
	#print numberoflines
	#print numberofobjects
	for q in range(numberofobjects):
		data.append([])
	#print data
	for i in range(numberoflines):
	        line = infile.readline()
	        if not line:
	            break
	        numbers = line.split()
		#print numbers
		for k in range(numberofobjects):
	        	data[k].append(float(numbers[k]))
			#print numbers[k]
	infile.close()
	#print data
    
	return data

def MagicNumbers(name):
	infile = open('%s.dat' % name, 'r')
	i = 0
	while True:
		line = infile.readline()
		if not line:
			break
		
		numbers = line.split()
		if not numbers:
			break
		q = len(numbers)
		#print numbers
		i = i + 1
	#print i, q
	return [i, q]

def Analytical(x, t):
	precision = 1000
	analytical = zeros(len(x))
	for i in range(len(analytical)):
		AB = 0
		A = 0
		B = 0
		for n in range(1, precision):
			A = -2*(-sin(pi*n) +pi*n)/(pi*pi*n*n)
			#for k in range(precision): 
			B = exp(-(n*n*pi*pi*t))#((-(n*n*pi*pi*t) )**k)/factorial(k)
			AB += A*B*sin(n*pi*x[i])
		analytical[i] = AB + (1-x[i])
	return analytical



l0 = 0.05
dt = (l0*l0)/(2)

name = 'data_time'
data = readfile(name)
#print data

anal_time = (len(data)/12)*dt
x = linspace(0, 1, MagicNumbers(name)[0])#readfile(name)[1]
print x
print len(data)
#print len(x)
#print data[:-1]
#print data[:]

# analytical at t
#plt.plot(x, Analytical(x, anal_time), '-r', label = 'analytical')
#plt.title('tidsutvikling')
#plt.hold('on')
##plt.plot(x, Analytical(x, 0.1), '-b')

#print data[-1]
plt.plot(x, data[int(len(data)/12)], '-b', label = 'simulation')
#plt.title('tidsutvikling')
plt.hold('on')
plt.plot(x, data[-1], '-r')

#plot all!
"""
for t in range(len(data)-1): # / number of dimensions in .dat file
	plt.plot(x, data[t], '-b', label = 'timestep%s' % t )
	plt.title('tidsutvikling')
	plt.hold('on')
"""


#plt.savefig('plot_%s.jpg' %name)
##plt.savefig('stable_crank_dt_0,04.jpg')

#plt.legend()
plt.show()

    
