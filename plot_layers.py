#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#plt.ion()
plt.rc('font', family='serif')

pol2 = np.loadtxt('run_layers_pol2.dat')
nt = np.loadtxt('run_layers_nt.dat')
nt=nt.transpose(); pol2=pol2.transpose()

plt.plot(nt[0],nt[1],'-or',label='NT-optimization')
plt.plot(pol2[0],pol2[1],'-ob',label='Pol2-optimization')
plt.xlim(0.8,120)
plt.semilogx()
plt.grid()
plt.xlabel('Number of layers')
plt.ylabel(r'$\sigma_t(\theta)$, mrad')
plt.legend()
plt.savefig('sigtang_vs_layers_varopt.pdf')
plt.show()

