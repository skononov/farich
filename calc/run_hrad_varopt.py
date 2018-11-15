#!/usr/bin/env python

import os, ROOT
from ROOT import gStyle, TFile, TH1F, TLegend, TCanvas


if __name__ == "__main__":
   os.system("./farichres -q -e 0.82 -p 3 -D 200 -T 40 -n 1.05 -N 4 -o mla4_fastopt.root data/mppc14160.dat")
   os.system("./farichres -q -e 0.82 -p 3 -D 200 -T 40 -n 1.05 -N 4 -mnt -o mla4_ntopt.root data/mppc14160.dat")
   os.system("./farichres -q -e 0.82 -p 3 -D 200 -T 40 -n 1.05 -N 4 -mpol2 -o mla4_pol2opt.root data/mppc14160.dat")
         
   f1 = TFile.Open("mla4_fastopt.root","READ")
   f2 = TFile.Open("mla4_ntopt.root","READ")
   f3 = TFile.Open("mla4_pol2opt.root","READ")

   gStyle.SetOptStat(0)
   gStyle.SetOptTitle(0)
   hrad1 = f1.Get("hrad")
   hrad1.SetLineWidth(2)
   hrad2 = f2.Get("hrad")
   hrad2.SetLineColor(2)
   hrad2.SetLineStyle(7)
   hrad2.SetLineWidth(2)
   hrad3 = f3.Get("hrad")
   hrad3.SetLineColor(4)
   hrad3.SetLineStyle(5)
   hrad3.SetLineWidth(2)
   
   c1 = TCanvas("c1","c1")
   
   l = TLegend(0.3,0.25,0.5,0.4)
   l.AddEntry(hrad1,"Fast","l")
   l.AddEntry(hrad2,"NT","l")
   l.AddEntry(hrad3,"Pol2","l")
   
   hrad1.Draw("l")
   hrad2.Draw("l same")
   hrad3.Draw("l same")
   l.Draw()
   
   c1.Print("hrad_mla4_varopt.pdf")

