#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TRandom3.h"
#include <gsl/gsl_sf_expint.h>

double func_sigmean(double* x, double *p)
{
    return sqrt((gsl_sf_expint_Ei(*x)-0.5772156649-log(*x))/(exp(*x)-1));
}

void sigmeantoymc(int nstat=10000) 
{
    TRandom3 rnd;
    
    int np = 11;
    double nmin = 0.1, nmax = 50., f = pow(nmax/nmin,1./(np-1));
    TGraph *gsigmean = new TGraph(np);
    
    TF1* fsigmean = new TF1("fsigmean", func_sigmean, 0., 1.1*nmax, 0);
    fsigmean->SetLineStyle(1);
    fsigmean->SetLineColor(kBlue);
    TF1* fsigmean0 = new TF1("fsigmean0", "1./sqrt(x)", 0., 1.1*nmax);
    fsigmean0->SetLineStyle(7);
    fsigmean0->SetLineColor(kRed);
    TF1* fsigmean1 = new TF1("fsigmean1", "1./sqrt(x/(1-exp(-x)))", 0., 1.1*nmax);
    fsigmean1->SetLineStyle(5);
    fsigmean1->SetLineColor(kMagenta);
    
    double n=nmin;
    for(int p=0; p<np; p++, n*=f) {
        double xmm = 0, x2mm = 0;
        int nmeas = 0;
        for(int i=0; i<nstat; i++) {
            int N = rnd.Poisson(n);
            if( N>0 ) {
                nmeas++;
                double xmean = 0;
                for(int j=0; j<N; j++)
                    xmean += rnd.Gaus(0,1)/N;
                xmm += xmean;
                x2mm += xmean*xmean;
            }
        }
        if( nmeas>0 ) {
            xmm /= nmeas;
            x2mm /= nmeas;
            gsigmean->SetPoint(p, n, sqrt(x2mm-xmm*xmm));
        }
    }
    
    TCanvas *c1 = new TCanvas("c1","c1");
    
    c1->DrawFrame(nmin/f,0.,nmax*f,1.05,";#nu;#sigma(#bar{x})");
    c1->SetLogx();
    gsigmean->Draw("p");
    fsigmean->Draw("same");
    fsigmean0->Draw("same");
    fsigmean1->Draw("same");
    
    TLegend *leg = new TLegend(0.15,0.15,0.6,0.4,"","brNDC");
    leg->SetTextFont(12);
    leg->SetTextSize(0.035);
    leg->AddEntry(gsigmean,"Monte Carlo","p");
    leg->AddEntry(fsigmean0,"1/#sqrt{#nu}","l");
    leg->AddEntry(fsigmean1,"1/#sqrt{#nu/(1-e^{-#nu})}","l");
    leg->AddEntry(fsigmean,"#sqrt{(Ei(#nu)-#gamma-log#nu)/(e^{#nu}-1)}","l");
    leg->Draw();
    
    c1->Print("sigmean.pdf");
}
