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
		for k in range(numberofobjects):
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
	

name = 'position_earth_euler'
data = readfile(name)
#print data
#x = linspace(0, ro_max, len(data[0]))#readfile(name)[1]
#print data[:-1]
#print data[:]

for t in range(len(data)/2): # / number of dimensions in .dat file
	plt.plot(data[0+t+t], data[1+t+t], '-b', label = 'object_%s' % t )
	plt.title('stability plot')
	plt.hold('on')

#plt.savefig('plot_%s.jpg' %name)
#plt.savefig('unstable_dt_0,09.jpg')

plt.show()
    
