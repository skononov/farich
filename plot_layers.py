#!/usr/bin/env python

import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

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
   
def clip_data(data, xlo=None, xup=None):
   mask = np.zeros_like(data, dtype=bool)
   if xlo is not None and xup is None:
      mask[0] = (data[0]<xlo)
   elif xlo is None and xup is not None:
      mask[0] = (data[0]>xup)
   elif xlo is not None and xup is not None:
      mask[0] = (data[0]<xlo or data[0]>xup)

   return ma.compress_cols(ma.array(data, mask=mask))
   
def main():
   if len(sys.argv) < 4:
      print "Usage: " + sys.argv[0] + " fast.dat nt.dat pol2.dat"
      return 1
      
   plt.rc('font', family='serif')


   fast = np.loadtxt(sys.argv[1],dtype='float32').transpose()
   nt = np.loadtxt(sys.argv[2],dtype='float32').transpose()
   pol2 = np.loadtxt(sys.argv[3],dtype='float32').transpose()

   nlmax = 30
   
   fast = clip_data(fast,xup=nlmax)
   nt = clip_data(nt,xup=nlmax)
   pol2 = clip_data(pol2,xup=nlmax)
   ymax = max(fast[1].max(),nt[1].max(),pol2[1].max())

   diff_nt_fast = diff(nt,fast,rel=True)
   diff_pol2_nt = diff(pol2,fast,rel=True)

   fig, (ax1, ax2) = plt.subplots(2,1)

   #ax1.set_title('Angle resolution for different methods')
   ax1.plot(fast[0],fast[1],'-ok',label='Fast')

   ax1.plot(nt[0],nt[1],'-sr',label='NT')

   ax1.plot(pol2[0],pol2[1],'-^b',label='Pol2')

   ax1.set_xlim(0.,nlmax*1.1)
   ax1.set_ylim(0.,ymax*1.1)
   #ax1.semilogx()
   ax1.grid()
   ax1.set_xlabel('Number of layers')
   ax1.set_ylabel(r'$\sigma_t(\theta)$, mrad')
   ax1.legend()
   ax1.text(0,ymax*1.13,'(a)')

   #ax2.set_title('Relative difference between methods')
   ax2.plot(diff_nt_fast[0], 100.*diff_nt_fast[1], '-o', label='NT-Fast')
   nd = (diff_pol2_nt[0]<=nlmax).sum()
   ax2.plot(diff_pol2_nt[0,:nd], 100.*diff_pol2_nt[1,:nd], '-s', label='Pol2-Fast')
   ax2.set_xlim(0.,nlmax*1.1)
   ax2.grid()
   ax2.set_xlabel('Number of layers')
   ax2.set_ylabel(r'Relative difference in $\sigma_t(\theta)$, %')
   ax2.legend()
   ylo, yup = ax2.get_ylim()
   ax2.text(0,yup+(yup-ylo)*0.06,'(b)')

   plt.tight_layout()

   fig.savefig('sigtang_vs_layers_varopt.pdf')
   #plt.show()

if __name__ == "__main__":
   rc = main()
   sys.exit(rc)


