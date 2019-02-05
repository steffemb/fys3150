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

name = 'file_10'
data = readfile(name)[0]
x = readfile(name)[1]

exact = 1.-(1-exp(-10))*x-exp(-10*x) 

plt.plot(x, data, '-b', label = 'my algoritm')
plt.title(name)
plt.hold('on')
plt.plot(x, exact, '-r', label = 'exact')
plt.legend()
plt.ylabel('u(x)')
plt.xlabel('x')
plt.savefig('plot_N_%s.jpg' %name)
plt.show()
    
