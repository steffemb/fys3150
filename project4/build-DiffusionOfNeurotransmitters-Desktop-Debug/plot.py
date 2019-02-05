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
		for k in range(numberofobjects-1):
	        	data[k].append(float(numbers[k]))
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
		#print numbers
		i = i + 1
	#print i, len(numbers)
	return [i, len(numbers)]
	

name = 'data_time'
data = readfile(name)
#print data
x = linspace(0, 1, len(data[0]))#readfile(name)[1]
#print data
#print len(data[0])
#print len(x)
#print data[:-1]
#print data[:]

plt.plot(x, data[len(data)/12], '-b')
plt.title('tidsutvikling')
plt.hold('on')
plt.plot(x, data[len(data)-4], '-r')
"""
for t in range(len(data)-1): # / number of dimensions in .dat file
	plt.plot(x, data[t], '-b', label = 'timestep%s' % t )
	plt.title('tidsutvikling')
	plt.hold('on')
"""
#plt.savefig('plot_%s.jpg' %name)
plt.savefig('stable_crank_dt_0,04.jpg')

plt.show()
    
