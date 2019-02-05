import matplotlib.pyplot as plt
from math import *
from numpy import*

def readfile(name):
    infile = open('%s.dat' % name, 'r')
    data = []
    x_list = []

    while True:
        line = infile.readline()
        if not line:
            break
        numbers = line.split()
        dat = float(numbers[0])
        x = float(numbers[1])
        #flux = float(numbers[2])
    
        data.append(dat)
        x_list.append(x)
        #fluxs.append(flux)
    infile.close()
    
    return [asarray(data), asarray(x_list)]

name = 'file_1000'
data = readfile(name)[0]
x = readfile(name)[1]
exact = 1.-(1-exp(-10))*x-exp(-10*x) 

error = zeros(len(data))
for i in range(len(data)):
    if data[i] == 0:
        error[i] = 0.
    else:
        error[i] = log10(abs( (exact[i]-data[i])/data[i] ))

print max(abs(error))

max_error = [-1.15000920684, -3.08918388604, -5.53390627447, -10.9243271902]
error_n = [10, 100, 1000, 10000]

plt.plot(error_n, max_error, '-b', label = 'error')
plt.title('error(n)')
plt.legend()
plt.ylabel('max relative error')
plt.xlabel('n')
plt.savefig('max_error.jpg')
plt.show()

"""
plt.plot(x, error, '-b', label = 'error')
plt.title(name + ' error')
plt.legend()
plt.ylabel('relative error')
plt.xlabel('x')
plt.savefig('error_n_1000.jpg')
plt.show()
"""
