#include <iostream>
#include <cassert>

#include "TCanvas.h"
#include "TStyle.h"
#include "TText.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"

TH1F* hpdf = nullptr;

double true_pdf(double *x, double *p)
{
    if (!hpdf) return 0.;

    return hpdf->Interpolate(x[0]-p[0]);
}

double erf_pdf(double *x, double *p)
{
    double r = (*x);

    return (TMath::Erf((r-p[0]+p[1]/2)/p[2]) - TMath::Erf((r-p[0]-p[1]/2)/p[3]))*p[4]*(1+p[5]*(r-p[0])+p[6]*(r-p[0])*(r-p[0]));
}

void unbinnedNLL_FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag) 
{
    npar = 7;

    if (flag==0) {
        f = ;
    }
    return 0;
}

void toysimrad(const char* fn, UInt_t nevents=10000)
{
    TFile* f = TFile::Open(fn, "READ");
    if (!f) return;

    TH1F* hrad = (TH1F*)f->Get("hrad");
    if (!hrad) return;

    hpdf = (TH1F*)hrad->Clone("hpdf");
    double binWidth = hrad->GetBinWidth(1);
    double npe = hrad->Integral()*binWidth;
    double rMean = hrad->GetMean();
    double rStdDev = hrad->GetStdDev();
    double rMin = hrad->GetXaxis()->GetXmin(), rMax = hrad->GetXaxis()->GetXmax();

    hpdf->Scale(1./npe);

    std::cout << "Toy simulation of FARICH for photoelectron density distribution given by " << fn << "\n"
              << "Mean radius = " << rMean << ", Radius std. dev. = " << rStdDev << ", Npe = " << npe << std::endl;

    TH1F* hradevent = (TH1F*)hrad->Clone("hradevent");

    TH1F* hradtrack1 = (TH1F*)hrad->Clone("hradtrack1");
    hradtrack1->Reset();
    hradtrack1->SetTitle("Mean radius per ring;ring radius, mm;number of events");

    TH1F* hradtrack2 = (TH1F*)hrad->Clone("hradtrack2");
    hradtrack2->Reset();
    hradtrack2->SetTitle("Ring radius fitted with max LH with true-like PDF;ring radius, mm;number of events");

    TH1F* hradtrack3 = (TH1F*)hrad->Clone("hradtrack3");
    hradtrack3->Reset();
    hradtrack3->SetTitle("Ring radius fitted with max LH with gaussian PDF;ring radius, mm;number of events");

    TF1* fTruePdf = new TF1("fTruePdf", true_pdf, 0., 2*rMean, 1);
    fTruePdf->SetLineColor(kBlue);
    fTruePdf->SetTitle("True-like PDF on hit radius;radius, mm");
    fTruePdf->SetParameter(0, 0.);

    TF1* fGausPdf = new TF1("fGausPdf", "gaus", 0., 2*rMean);
    fGausPdf->SetLineColor(kRed);
    fGausPdf->SetTitle("Gaussian PDF on hit radius;radius, mm");
    fGausPdf->SetParameters(1/sqrt(2*M_PI)/rStdDev, rMean, rStdDev);
    fGausPdf->FixParameter(0, fGausPdf->GetParameter(0));
    fGausPdf->FixParameter(2, fGausPdf->GetParameter(2));
    fGausPdf->SetNpx(1000);

    TF1* fErfPdf = new TF1("fErfPdf", erf_pdf, 0., 2*rMean, 7);
    fErfPdf->SetParameters(rMean, 4*rStdDev, 1, 1, 0.1, -0.1, 0.01);
    fErfPdf->SetLineColor(kGreen+2);
    fErfPdf->SetNpx(1000);
    //fit shape to the real pdf
    Int_t fitStatus = hpdf->Fit(fErfPdf, "N"); 
    if (fitStatus!=0) {
        std::cerr << "Error fitting erf PDF to the real PDF" << std::endl;
        return;
    }
    //fix shape and leave X position free
    for(int ipar = 1; ipar < fErfPdf->GetNpar(); ipar++)
        fErfPdf->FixParameter(ipar, fErfPdf->GetParameter(ipar));
    //fErfPdf->ReleaseParameter(4);
    double meanErfPdfShift = fErfPdf->Mean(rMin, rMax) - fErfPdf->GetParameter(0);
    double erfPdfNormCoef = fErfPdf->GetParameter(4)/fErfPdf->Integral(rMin, rMax);

    TCanvas *c1 = new TCanvas("c1","c1",1002,1028);
    c1->Divide(2,2);

    c1->cd(1);
    hpdf->Draw("hist");
    fTruePdf->SetNpx(1000);
    fTruePdf->Draw("same");
    fGausPdf->Draw("same");
    fErfPdf->Draw("same");

    TText *text = new TText(0.5,0.35,TString::Format("Npe=%.1f",npe));
    text->SetNDC();
    text->Draw();
    c1->Update();

    //return;

    c1->cd(2);
    TRandom3 grnd(1234);

    for(UInt_t i=0; i<nevents; i++) {
        Int_t n = grnd.Poisson(npe);
        hradevent->Reset();
        
        fTruePdf->SetParameter(0, 0.);
        hradevent->FillRandom("fTruePdf", n);

        hradtrack1->Fill(hradevent->GetMean());
/*      
        fTruePdf->SetParameter(0, hradevent->GetMean()-rMean);
        fitStatus = hradevent->Fit(fTruePdf, "LNQ", "goff");
        if (fitStatus==0) 
            hradtrack2->Fill(rMean+fTruePdf->GetParameter(0));
*/
        fErfPdf->SetParameter(0, hradevent->GetMean()-meanErfPdfShift);
        //fErfPdf->SetParameter(4, n*erfPdfNormCoef*binWidth);
        fitStatus = hradevent->Fit(fErfPdf, "LNQ", "goff");
        if (fitStatus==0) 
            hradtrack2->Fill(fErfPdf->GetParameter(0)+meanErfPdfShift);

        fGausPdf->SetParameter(1, hradevent->GetMean());
        //fGausPdf->SetParameter(0, n/sqrt(2*M_PI)/rStdDev*binWidth);
        fitStatus = hradevent->Fit(fGausPdf, "LNQB+", "goff");
        if (fitStatus==0)
            hradtrack3->Fill(fGausPdf->GetParameter(1));

        //hradevent->Draw();
        //gPad->Update();

        if(i%1000==0) std::cout << "Event " << i << std::endl;
    }

    fTruePdf->SetParameter(0, 0.);

    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(10111);

    c1->cd(2);
    double rmin = hradtrack1->GetMean()-5*hradtrack1->GetRMS(), rmax = hradtrack1->GetMean()+5*hradtrack1->GetRMS();
    hradtrack1->SetAxisRange(rmin, rmax);
    hradtrack1->Fit("gaus","I","",rmin,rmax);

    c1->cd(3);
    rmin = hradtrack2->GetMean()-5*hradtrack2->GetRMS(); rmax = hradtrack2->GetMean()+5*hradtrack2->GetRMS();
    hradtrack2->SetAxisRange(rmin, rmax);
    hradtrack2->Fit("gaus","I","",rmin,rmax);

    c1->cd(4);
    rmin = hradtrack3->GetMean()-5*hradtrack3->GetRMS(); rmax = hradtrack3->GetMean()+5*hradtrack3->GetRMS();
    hradtrack3->SetAxisRange(rmin, rmax);
    hradtrack3->Fit("gaus","I","",rmin,rmax);
}