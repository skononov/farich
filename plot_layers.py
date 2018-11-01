#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

#plt.ion()
plt.rc('font', family='serif')

def diff(a, b, rel=False):
   d = []
   ita = np.nditer(a[0],flags=['multi_index'])
   itb = np.nditer(b[0],flags=['multi_index'])
   while not (ita.finished or itb.finished):
      na = ita[0]
      nb = itb[0]
      if na == nb:
         i = ita.multi_index[0]
         j = itb.multi_index[0]
         if rel:
            dab = (a[1,i]-b[1,j])/a[1,i]
         else:
            dab = a[1,i]-b[1,j]
         d.append((na,dab))
         ita.iternext()
         itb.iternext()
      elif na > nb:
         itb.iternext()
      else:
         ita.iternext()
   return np.array(d,dtype='float').transpose()

fast = np.loadtxt('run_layers_fast.dat')
nt = np.loadtxt('run_layers_nt.dat')
pol2 = np.loadtxt('run_layers_pol2.dat')
fast = fast.transpose()
nt = nt.transpose()
pol2 = pol2.transpose()

diff_nt_fast = diff(nt,fast,rel=True)
diff_pol2_nt = diff(pol2,fast,rel=True)

fig, (ax1, ax2) = plt.subplots(2,1)

#ax1.set_title('Angle resolution for different methods')
ax1.plot(fast[0],fast[1],'-ok',label='Fast')
ax1.plot(nt[0],nt[1],'-sr',label='NT')
ax1.plot(pol2[0][:9],pol2[1][:9],'-^b',label='Pol2')
ax1.set_xlim(0.,22.)
ax1.set_ylim(0.,2.2)
#ax1.semilogx()
ax1.grid()
ax1.set_xlabel('Number of layers')
ax1.set_ylabel(r'$\sigma_t(\theta)$, mrad')
ax1.legend()
ax1.text(18.,1.8,'(a)')

#ax2.set_title('Relative difference between methods')
ax2.plot(diff_nt_fast[0], 100.*diff_nt_fast[1], '-o', label='NT-Fast')
ax2.plot(diff_pol2_nt[0], 100.*diff_pol2_nt[1], '-s', label='Pol2-Fast')
ax2.grid()
ax2.set_xlabel('Number of layers')
ax2.set_ylabel(r'Relative difference in $\sigma_t(\theta)$, %')
ax2.legend()
ax2.text(8.,.55,'(b)')

plt.tight_layout()

fig.savefig('sigtang_vs_layers_varopt.pdf')
plt.show()

