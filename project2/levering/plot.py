import matplotlib.pyplot as plt
from math import *
from numpy import*

def readfile(name):
    infile = open('%s.dat' % name, 'r')
    eigvec_1 = []
    eigvec_2 = []
    eigvec_3 = []
    eigvec_1.append(0)
    eigvec_2.append(0)
    eigvec_3.append(0)

    while True:
        line = infile.readline()
        if not line:
            break
        numbers = line.split()
        eigvec_1_ = float(numbers[0])
        eigvec_2_ = float(numbers[1])
        eigvec_3_ = float(numbers[2])
	ro_max = float(numbers[3])
    
        eigvec_1.append(eigvec_1_)
        eigvec_2.append(eigvec_2_)
        eigvec_3.append(eigvec_3_)
    infile.close()
    
    return [asarray(eigvec_1), asarray(eigvec_2), asarray(eigvec_3), ro_max]

name = 'eigenvec'
data = readfile(name)
ro_max = data[3]
x = linspace(0, ro_max, len(data[0]))#readfile(name)[1]
#print data[:-1]
#print data[:]

plt.plot(x, data[0], '-b', label = 'psi_1')
plt.title(name)
plt.hold('on')
plt.plot(x, data[1], '-r', label = 'psi_2')
plt.hold('on')
plt.plot(x, data[2], '-g', label = 'psi_3')
plt.legend()
plt.ylabel('psi(x)')
plt.xlabel('ro')
#plt.savefig('plot_%s.jpg' %name)
plt.savefig('plot_omega_5.jpg')
plt.show()
    
