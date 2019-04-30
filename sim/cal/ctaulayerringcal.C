#include "TCanvas.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../analysis/aerdis.C"
#include "../analysis/particles.h"
#include "getflist.h"

const char *CT_DEFAULT = "\033[0m";
const char *CT_INVERSE = "\033[7m";
const char *CT_BGREDBOLD = "\033[41;1m";
const char *CT_FGCYANBOLD = "\033[1;36m";
const char *CT_FGREDBOLD = "\033[1;31m";

const Int_t nLayers = 4;
const Int_t nBins = 45;
const Int_t nCosPars = 5;
const Int_t nTotPars = nCosPars * 3 * nLayers;
const Double_t deg2rad = TMath::Pi() / 180.;
const Double_t sqrt2pi = sqrt(2 * TMath::Pi());
const Double_t riMax = 1.07;
const Float_t probMin = 0.1;
Float_t probCutToDraw = 0.01;
Bool_t fitPDF = kFALSE;
Bool_t DEBUG = kFALSE;
Float_t pSelected = 0., thetaSelected = 0.;

Double_t func_pdf2d(Double_t *x, Double_t *par)
{
    Double_t phic = x[0];
    Double_t thetac = x[1];

    if (thetac < 0 || thetac > 1570 || phic < -180 || phic > 180) {
        TF1::RejectPoint();
        return 0.;
    }

    Double_t C = par[0], M = par[nCosPars], S = par[2 * nCosPars];
    Double_t cosnphic = 0;

    for (Int_t icp = 1; icp < nCosPars; icp++) {
        cosnphic = cos(icp * deg2rad * phic);
        C += par[icp] * cosnphic;
        M += par[nCosPars + icp] * cosnphic;
        S += par[2 * nCosPars + icp] * cosnphic;
    }

    if (C <= 0 || S <= 0)
        return 0.;

    return C / sqrt2pi / S * exp(-(thetac - M) * (thetac - M) / (2 * S * S));
}

int ctaulayerringcal(TString fnpat = "/home/skononov/geant4/farich/runs/ctau/reco/"
                                     "pi_ctau_mppc2_px3_d200_mla4_*deg_*gev.root",
                     TString outfn = "fitpar_layers.dat")
{
    Ssiz_t pos = fnpat.Last('/') + 1, len = fnpat.Index("_", pos) - pos;
    TString particle = fnpat(pos, len);

    Double_t M = get_particle_mass(particle) / 1e3;
    if (M == 0)
        return -1;

    TList *flist = getflist(fnpat);
    if (!flist)
        return -1;

    cout << flist->GetSize() << " input files defined by the pattern" << endl;

    TIter piter(flist);
    ParFileRef *ref;
    TH1F *hang, *hnpe;
    TTree *th;

    if (pSelected > 0) {
        while ((ref = (ParFileRef *)piter())) {
            if (ref->GetParameter(0) != thetaSelected || ref->GetParameter(1) != pSelected)
                flist->Remove(ref);
        }
        if (flist->GetSize() == 0) {
            cout << "!No files for selected p=" << pSelected << " and theta=" << thetaSelected << endl;
            return -1;
        }
        piter.Reset();
    }

    // Determine lower point on momentum
    Double_t Plo = 0;
    while ((ref = (ParFileRef *)piter())) {
        Double_t Theta = ref->GetParameter(0);
        Double_t P = ref->GetParameter(1);
        if (DEBUG)
            cout << "Theta=" << Theta << " P=" << P << " " << ref->GetName() << endl;

        if (P <= Plo)
            continue;

        ref->GetObject("hang", hang);

        if (hang->GetEntries() < 2000)
            Plo = P;
    }
    piter.Reset();
    if (Plo == 0)
        cout << "No files rejected with statistics cut" << endl;
    else
        cout << "Reject files with P<=" << setprecision(5) << Plo << " GeV/c because of low statistics" << endl;

    if (fitPDF)
        cout << "Fit Cherenkov angle distributions" << endl;
    else
        cout << "Do not fit Cherenkov angle distributions" << endl;

    //	TF1* fgaus=new
    // TF1("fgaus","[0]/2.506628/[2]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",0,400);
    //	fgaus->SetParError(1,10);
    //	fgaus->SetParError(2,1);
    TF1 *fcos = new TF1("fcos",
                        "[0]+[1]*cos(0.01745329*x)+[2]*cos(2*0.01745329*x)+[3]*"
                        "cos(3*0.01745329*x)+[4]*cos(4*0.01745329*x)",
                        -180, 180);
    TF2 *fpdf2d = new TF2("fpdf2d", func_pdf2d, -180, 180, 0, 400, 3 * nCosPars);

    const char *good = " (fit success)";
    const char *failed = " \033[41;1m(fit failed)\033[0m";

    Float_t pars[nTotPars];
    const TString expr = TString::Format("angle:alpha>>haa(%d,-180,180,%d,0,400)", nBins, fitPDF ? 400 : 800);
    const char *names[4] = {"Constant", "Mean", "Sigma", "Chisquare"};
    const Int_t kNotDraw = 1 << 9;
    TObjArray fitarr;
    TGraphErrors *gpar[4][nLayers];
    memset(gpar, 0, sizeof(gpar));
    Bool_t draw[4][nLayers];
    Double_t prevTheta = -1;

    gStyle->SetOptStat(10);
    gStyle->SetOptFit(1011);
    TCanvas *c1 = new TCanvas("c1", "c1", 260, 500, 500, 440);
    TCanvas *c2 = new TCanvas("c2", "c2", 260, 50, 500, 440);
    TCanvas *c3 = new TCanvas("c3", "c3", 780, 50, 500, 440);
    TCanvas *c4 = new TCanvas("c4", "c4", 780, 500, 500, 440);
    TCanvas *canvases[4] = {c1, c2, c3, c4};

    ofstream fout(outfn);
    fout.precision(6);

    while ((ref = (ParFileRef *)piter())) {
        ref->cd();
        Double_t Theta = ref->GetParameter(0);
        Double_t P = ref->GetParameter(1);
        if (P <= Plo)
            continue;

        cout.precision(4);
        cout << "Dip angle=" << Theta << " deg, momentum=" << P << " GeV/c" << endl;

        Double_t Beta = 1 / sqrt(1 + M * M / P / P);

        ref->GetObject("th", th);
        if (!th) {
            cout << CT_BGREDBOLD << "No th tree found!" << CT_DEFAULT << endl;
            break;
        }

        ref->GetObject("hnpe", hnpe);
        Double_t nevents = hnpe->GetEntries();

        memset(pars, 0, sizeof(pars));
        memset(draw, 0, sizeof(draw));

        for (Int_t ig = 0; ig < 3; ig++) {
            for (Int_t il = 0; il < nLayers; il++) {
                if (gpar[ig][il]) {
                    delete gpar[ig][il];
                    gpar[ig][il] = 0;
                }
            }
        }

        gDirectory->Delete("pa*");

        for (Int_t il = 0; il < nLayers; il++) {
            gSystem->ProcessEvents();
            TString hname = TString::Format("hang%d", il + 1);
            ref->GetObject(hname, hang);
            if (!hang) {
                cout << CT_BGREDBOLD << "No " << hname << " histogram!" << CT_DEFAULT << endl;
                continue;
            }
            if (hang->GetEntries() == 0) {
                cout << " Skip other layers because of no hits" << endl;
                break;
            }

            th->Draw(expr, TString::Format("layer==%d", il + 1), "goff");
            TH2F *haa = 0;
            ref->GetObject("haa", haa);
            if (!haa) {
                cout << CT_BGREDBOLD << "No haa histogram!" << CT_DEFAULT << endl;
                break;
            }

            if (haa->GetEntries() < 10 * nBins) {
                if (haa->GetEntries() < 3 * nBins) {
                    cout << " Skip other layers because of low number of hits" << endl;
                    break;
                }
                haa->RebinX(3);
            }
            Double_t xbinw = haa->GetXaxis()->GetBinWidth(1);
            Double_t ybinw = haa->GetYaxis()->GetBinWidth(1);

            Axis_t ymin = haa->GetYaxis()->GetBinLowEdge(haa->FindFirstBinAbove(0, 2));
            Axis_t ymax = haa->GetYaxis()->GetBinUpEdge(haa->FindLastBinAbove(0, 2));

            TString pname("pa");
            pname += il + 1;
            TProfile *pa = haa->ProfileX(pname, -1, -1, "s");
            Int_t np = 0;
            for (Int_t ig = 0; ig < 3; ig++)
                gpar[ig][il] = new TGraphErrors(nBins);
            for (Int_t bin = 1; bin <= nBins; bin++) {
                Double_t ent = pa->GetBinEntries(bin);
                if (ent < 5)
                    continue;
                Double_t sqrtent = sqrt(ent);
                Double_t x = pa->GetBinCenter(bin);
                Double_t val = pa->GetBinContent(bin);
                Double_t rms = pa->GetBinError(bin);
                gpar[0][il]->SetPoint(np, x, ent / xbinw / nevents);
                gpar[0][il]->SetPointError(np, 0, sqrtent / xbinw / nevents);
                gpar[1][il]->SetPoint(np, x, val);
                gpar[1][il]->SetPointError(np, 0, rms / sqrtent);
                gpar[2][il]->SetPoint(np, x, rms);
                gpar[2][il]->SetPointError(np, 0, rms / sqrt(2.) / sqrtent);
                np++;
            }
            for (Int_t ig = 0; ig < 3; ig++) {
                gpar[ig][il]->Set(np);
                fcos->SetParameter(0, TMath::Mean(np, gpar[ig][il]->GetY()));
                if (Theta > 0)
                    fcos->SetParameter(1, 0.);
                else
                    fcos->FixParameter(1, 0.);
                for (Int_t icp = 2; icp < nCosPars; icp++) {
                    fcos->FixParameter(icp, 0.);
                }
                gpar[ig][il]->Fit(fcos, "QR0", "goff");
            }

            cout << " Layer " << il + 1 << ". ";
            Int_t rc = 0;
            if (fitPDF) {
                for (Int_t ig = 0; ig < 3; ig++) {
                    Double_t val = gpar[ig][il]->GetFunction("fcos")->GetParameter(0);
                    if (ig == 0) // special value for gaus Const parameter
                        val *= xbinw * nevents;
                    fpdf2d->SetParameter(ig * nCosPars, val);
                    fpdf2d->SetParError(ig * nCosPars, 0.1 * val);
                    if (Theta > 0) {
                        val = gpar[ig][il]->GetFunction("fcos")->GetParameter(1);
                        if (ig == 0) // special value for gaus Const parameter
                            val *= xbinw * nevents;
                        fpdf2d->SetParameter(ig * nCosPars + 1, val);
                        fpdf2d->SetParError(ig * nCosPars + 1, 0.1 * val);
                    } else
                        fpdf2d->FixParameter(ig * nCosPars + 1, 0.);
                    for (Int_t icp = 2; icp < nCosPars; icp++) {
                        fpdf2d->FixParameter(ig * nCosPars + icp, 0.);
                    }
                }

                fpdf2d->SetRange(-180., ymin, 180., ymax);

                if (Theta > 0) {
                    for (Int_t icp = 1; icp < nCosPars; icp++) {
                        fpdf2d->ReleaseParameter(icp);
                        fpdf2d->ReleaseParameter(nCosPars + icp);
                        fpdf2d->ReleaseParameter(2 * nCosPars + icp);
                        rc = haa->Fit(fpdf2d, "LQR0", "goff");
                        if (fpdf2d->GetProb() > probMin)
                            break;
                    }
                } else {
                    rc = haa->Fit(fpdf2d, "LQR0", "goff");
                }

                cout << "pdf2d fit: prob=" << (fpdf2d->GetProb() < probCutToDraw ? CT_FGREDBOLD : "")
                     << fpdf2d->GetProb() << CT_DEFAULT << (rc ? failed : good) << endl;

                // Fill pars array and fcos functions with fitted parameters and mark
                // for drawing
                for (Int_t ig = 0; ig < 3; ig++) {
                    for (Int_t icp = 0; icp < nCosPars; icp++) {
                        Double_t val = fpdf2d->GetParameter(ig * nCosPars + icp);
                        if (ig == 0) // special treating for gaus Const parameter
                            val *= 1. / xbinw / ybinw / nevents;
                        gpar[ig][il]->GetFunction("fcos")->SetParameter(icp, val);
                        pars[nCosPars * (3 * il + ig) + icp] = val;
                    }

                    if (!DEBUG && fpdf2d->GetProb() > probCutToDraw)
                        draw[ig][il] = 0;
                    else {
                        draw[ig][il] = 1;
                        gpar[ig][il]->GetFunction("fcos")->ResetBit(kNotDraw);
                        gpar[ig][il]->GetFunction("fcos")->SetLineColor(kRed + 3 - il);
                        gpar[ig][il]->SetLineColor(kBlue + 3 - il);
                        gpar[ig][il]->SetMarkerColor(kBlue + 3 - il);
                    }
                }

#if 0
//				haa->FitSlicesY((TF1*)0,1,nBins,5,"QNRL",&fitarr);
				fgaus->SetRange(ymin,ymax);
				fgaus->SetParameters(haa->GetEntries()/nBins,(ymin+ymax)/2,ymax-ymin);
				fgaus->SetParError(0,0.1*fgaus->GetParameter(0));

				haa->FitSlicesY(fgaus,1,nBins,5,"QNRL",&fitarr);

				TH1D *hchi2=(TH1D*)fitarr[3];

				for(Int_t ig=0; ig<3; ig++) {
					TH1D *hpar=(TH1D*)fitarr[ig];
					if( !hpar ) {
	                    cout<<CT_BGREDBOLD<<"No slice fit histogram for "<<names[ig]<<CT_DEFAULT<<endl;
						break;
	                }

					Int_t np=0;
					for(Int_t bin=1; bin<=hpar->GetNbinsX(); bin++) {
						Double_t val=hpar->GetBinContent(bin), err=hpar->GetBinError(bin);
						if( val<0 || err<=0 ) continue;
						if( hchi2->GetBinContent(bin)>1.7 ) continue;
						if( ig==0 ) {
							val*=1./xbinw/ybinw/nevents;
							err*=1./xbinw/ybinw/nevents;
                    	}
						gpar[ig][il]->SetPoint(np,hpar->GetBinCenter(bin),val);
						gpar[ig][il]->SetPointError(np,0,err);
						np++;
					}
					gpar[ig][il]->Set(np);
				} //for(Int_t ig=0; ig<3; ig++)
				gpar[3][il] = new TGraphErrors(hchi2->GetNbinsX());
				gpar[3][il]->SetMarkerColor(kBlue+3-il);
				for(Int_t bin=1; bin<=hchi2->GetNbinsX(); bin++) {
					gpar[3][il]->SetPoint(bin-1,hchi2->GetBinCenter(bin),hchi2->GetBinContent(bin));
				}
				fitarr.Clear();
#endif
            } else { // if( !fitPDF )
                for (Int_t ig = 0; ig < 3; ig++) {
                    if (Theta > 0) {
                        for (Int_t icp = 1; icp < nCosPars; icp++) {
                            fcos->ReleaseParameter(icp);
                            rc = gpar[ig][il]->Fit(fcos, "Q0R", "goff");
                            if (fcos->GetProb() > probMin)
                                break;
                        }
                    } else {
                        rc = gpar[ig][il]->Fit(fcos, "Q0R", "goff");
                    }
                    cout << " " << names[ig] << " fit: prob=" << (fcos->GetProb() < 0.01 ? CT_FGREDBOLD : "")
                         << fcos->GetProb() << CT_DEFAULT << (rc ? failed : good);
                    for (Int_t icp = 0; icp < nCosPars; icp++)
                        pars[nCosPars * (3 * il + ig) + icp] = fcos->GetParameter(icp);

                    if (!DEBUG && fcos->GetProb() > 0.01)
                        draw[ig][il] = 0;
                    else {
                        draw[ig][il] = 1;
                        draw[3][il] = 1;
                        gpar[ig][il]->GetFunction("fcos")->ResetBit(kNotDraw);
                        gpar[ig][il]->GetFunction("fcos")->SetLineColor(kRed + 3 - il);
                        gpar[ig][il]->SetLineColor(kBlue + 3 - il);
                        gpar[ig][il]->SetMarkerColor(kBlue + 3 - il);
                    }
                } // for(Int_t ig=0; ig<3; ig++)
                cout << endl;
            } // if( fitPDF )
        }     // for(Int_t il=0; il<nLayers; il++)

        // Draw TGraphs along with fitted functions
        if (!gpar[0][0] || !gpar[1][0] || !gpar[2][0]) {
            cout << CT_FGREDBOLD << "No parameter histogram for the first layer" << CT_DEFAULT << endl;
        } else {
            for (Int_t ig = 0; ig < 3; ig++) {
                Double_t ymin = 0, ymax = 0;
                for (Int_t il = 0; il < nLayers; il++) {
                    if (!gpar[ig][il])
                        break;
                    if (!draw[ig][il])
                        continue;
                    if (ymin == 0 && ymax == 0)
                        ymin = ymax = gpar[ig][il]->GetY()[0];
                    for (Int_t i = 0; i < gpar[ig][il]->GetN(); i++) {
                        if (ymin > gpar[ig][il]->GetY()[i] - gpar[ig][il]->GetEY()[i])
                            ymin = gpar[ig][il]->GetY()[i] - gpar[ig][il]->GetEY()[i];
                        if (ymax < gpar[ig][il]->GetY()[i] + gpar[ig][il]->GetEY()[i])
                            ymax = gpar[ig][il]->GetY()[i] + gpar[ig][il]->GetEY()[i];
                    }
                }
                if (ymin < ymax) {
                    canvases[ig]->cd();
                    canvases[ig]->DrawFrame(-190, ymin - 0.1 * (ymax - ymin), 190, ymax + 0.1 * (ymax - ymin),
                                            TString::Format("%s: #theta=%d^{0} P=%6.4gGeV/c;#alpha, deg;%s", names[ig],
                                                            (Int_t)Theta, P, names[ig]));
                    for (Int_t il = 0; il < nLayers; il++) {
                        if (!gpar[ig][il])
                            break;
                        if (!draw[ig][il])
                            continue;
                        gpar[ig][il]->Draw("pe");
                    }
                    canvases[ig]->Update();
                }
            }
        }

        cout.precision(4);

        if (Theta != prevTheta) { // first file with next Theta
            // Put a line with zero parameters at threshold beta
            cout << " " << 1 / riMax << " " << Theta << "  \033[1;33m zero line at threshold" << CT_DEFAULT << endl;
            fout << 1 / riMax << " " << Theta;
            for (Int_t ipar = 0; ipar < nTotPars; ipar++) {
                fout << " 0";
            }
            fout << endl;
            prevTheta = Theta;
        }

        cout << " " << Beta << " " << Theta;
        fout << Beta << " " << Theta;
        for (Int_t ipar = 0; ipar < nTotPars; ipar++) {
            if (ipar % nCosPars == 0)
                cout << CT_DEFAULT;
            if (ipar % (2 * nCosPars) == 0)
                cout << "  \033[1;33m";
            cout << " " << pars[ipar];
            fout << " " << pars[ipar];
        }
        cout << CT_DEFAULT << endl;
        fout << endl;

        ref->Close();
    }

    fout.close();

    flist->Delete();
    delete flist;
}
