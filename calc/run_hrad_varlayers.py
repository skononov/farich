#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, re, ROOT
import  subprocess as sp
from ROOT import gStyle, TFile, TH1F, TLegend, TCanvas

if __name__ == "__main__":
   nlarr = (1, 2, 4, 10, 20)
   styles = [1, 2, 7, 5, 9]
   #cols = [1, 2, 4, 6, 9]
   cols = [1, 1, 1, 1, 1]
   
   fnames = ["mla%d_d200_t40.root" % nl for nl in nlarr]
   
   sigrad1 = []
   sigradt = []
   sigang1 = []
   sigangt = []
   
   resigrad1 = re.compile(r'Ошибка радиуса на 1 фотон.*: +([0-9.]+) мм +\(([0-9.]+) мм\)')
   resigradt = re.compile(r'Ошибка радиуса на трек.*: +([0-9.]+) мм +\(([0-9.]+) мм\)')
   resigang1 = re.compile(r'Ошибка угла на 1 фотон.*: +([0-9.]+) мрад +\(([0-9.]+) мрад\)')
   resigangt = re.compile(r'Ошибка угла на трек.*: +([0-9.]+) мрад +\(([0-9.]+) мрад\)')
   
   def parse_sigma(s,resig):
      m = resig.search(s)
      if m is None or len(m.groups())<2:
         print 'Error: Can not find RE ' + resig
         return None
      return (float(m.group(1)),float(m.group(2)))
         

   for nl,fn in zip(nlarr,fnames):
      print '\nNumber of layers %d' % nl
      cmd = "./farichres -q -e 0.82 -p 3 -D 200 -T 40 -n 1.05 -N %d -mpol2 -o %s data/mppc14160.dat" % (nl,fn)
      print 'Running: ' + cmd
      output = sp.check_output(cmd,stderr=sp.STDOUT,shell=True)

      sigrad1.append(parse_sigma(output,resigrad1))
      sigradt.append(parse_sigma(output,resigradt))
      sigang1.append(parse_sigma(output,resigang1))
      sigangt.append(parse_sigma(output,resigangt))

   sigrad1 = zip(*sigrad1)
   sigradt = zip(*sigradt)
   sigang1 = zip(*sigang1)
   sigangt = zip(*sigangt)
   
   files = [TFile.Open(fn,"READ") for fn in fnames] 
   hists = [f.Get("hrad") for f in files]

   gStyle.SetOptStat(0)
   gStyle.SetOptTitle(0)
   c1 = TCanvas("c1","c1")

   cmax = max([h.GetMaximum() for h in hists])
   xmin = min([h.GetXaxis().GetXmin() for h in hists])
   xmax = max([h.GetXaxis().GetXmax() for h in hists])
   c1.DrawFrame(xmin-0.05*(xmax-xmin),0.,xmax+0.05*(xmax-xmin),cmax*1.1,\
      h.GetTitle()+';'+h.GetXaxis().GetTitle()+';'+h.GetYaxis().GetTitle())
   l = TLegend(0.3,0.5,0.85,0.85)
   l.SetBorderSize(0)
   l.SetTextFont(42)
   
   radarr = []
   npearr = []

   for i,nl in enumerate(nlarr):
      hists[i].SetLineWidth(2)
      hists[i].SetLineStyle(styles[i])
      hists[i].SetLineColor(cols[i])
   
      radarr.append(hists[i].GetMean())
      npearr.append(hists[i].Integral()*hists[i].GetBinWidth(1))
      
      l.AddEntry(hists[i],"N=%d, #sigma_{1}(R)=%4.2gmm, #sigma_{t}(#theta)=%4.2gmrad" % (nl,sigrad1[0][i],sigangt[0][i]),"l")
      
      hists[i].Draw("l same")

   l.Draw()
   
   c1.Print("hrad_varlayers_pol2opt.pdf")

   print r'Число слоев & ' + ' & '.join(['%2d' % nl for nl in nlarr]) + r'\\\hline'
   print r'$N_\mathrm{pe}$ & ' + ' & '.join(['%4.2g' % npe for npe in npearr]) + r'\\\hline'
   print r'$R$, мм & ' + ' & '.join(['%5.3g' % r for r in radarr]) + r'\\\hline'
   print r'\multicolumn{6}{c|}{Ошибки без учета размера пикселя}\\\hline'
   print r'$\sigma_1(R)$, мм & ' + ' & '.join(['%4.2g' % s for s in sigrad1[0]]) + r'\\\hline'
   print r'$\sigma_t(R)$, мм & ' + ' & '.join(['%4.2g' % s for s in sigradt[0]]) + r'\\\hline'
   print r'$\sigma_1(\theta)$, мрад & ' + ' & '.join(['%4.2g' % s for s in sigang1[0]]) + r'\\\hline'
   print r'$\sigma_t(\theta)$, мрад & ' + ' & '.join(['%4.2g' % s for s in sigangt[0]]) + r'\\\hline'
   print r'\multicolumn{6}{c|}{Ошибки с учетом размера пикселя}\\\hline'
   print r'$\sigma_1(R)$, мм & ' + ' & '.join(['%4.2g' % s for s in sigrad1[1]]) + r'\\\hline'
   print r'$\sigma_t(R)$, мм & ' + ' & '.join(['%4.2g' % s for s in sigradt[1]]) + r'\\\hline'
   print r'$\sigma_1(\theta)$, мрад & ' + ' & '.join(['%4.2g' % s for s in sigang1[1]]) + r'\\\hline'
   print r'$\sigma_t(\theta)$, мрад & ' + ' & '.join(['%4.2g' % s for s in sigangt[1]]) + r'\\\hline'

   try:
      input('Press Enter...')
   except:
      pass

