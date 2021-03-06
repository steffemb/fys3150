import matplotlib.pyplot as plt
from math import *
from numpy import*

my_algo = [0, 0, 0, 0]
LU_algo = [0, 0, 0.05, 5.84]
n = [10, 100, 1000, 10000]


plt.plot(n, my_algo, '-b', label = 'my algoritm')
plt.title('bencmark of algoritms')
plt.hold('on')
plt.plot(n, LU_algo, '-r', label = 'LU_algo')
plt.legend()
plt.ylabel('time [s]')
plt.xlabel('n')
plt.savefig('benchmark.jpg')
plt.show()
    
