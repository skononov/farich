{
#include <iomanip>

    const char *CT_DEFAULT = "\033[0m";
    const char *CT_INVERSE = "\033[7m";
    const char *CT_BGREDBOLD = "\033[41;1m";
    const char *CT_FGCYANBOLD = "\033[1;36m";
    const char *CT_FGREDBOLD = "\033[1;31m";

    if (!gInterpreter->IsLoaded("libMyRoot.so"))
        gSystem->Load("libMyRoot.so");

    TList *flist = getflist("/home/skononov/geant4/farich/runs/ctau/reco/"
                            "pi_ctau_mppc2_px3_d200_mla4_*deg_*gev.root");
    if (!flist)
        return;

    cout << flist->GetSize() << " files found" << endl;

    const Double_t M = 0.13957; // pion mass GeV/c
    const Double_t CL = 0.95;   // confidence level for rejection

    TIter piter(flist);
    ParFileRef *ref;
    TProfile *pang;

#if 0
	while( (ref=(ParFileRef*)piter()) )
	{
		if( ref->GetParameter(0)!=16. || ref->GetParameter(1)!=0.8543 )
			flist->Remove(ref);

	}
	if( flist->GetSize()==0 ) {
		cout<<" No files found"<<endl;
		return;
	}
    piter.Reset();
#endif

    ref = (ParFileRef *)flist->At(0);
    ref->GetObject("pangalp", pang);

    const Int_t nbins = pang->GetNbinsX();
    TGraphErrors *gang = new TGraphErrors(nbins);
    TGraphErrors *gsig = new TGraphErrors(nbins);

    const Double_t CL1 = pow(CL, 1. / nbins); // rejection confidence level for single point
    const Double_t devLimit = ROOT::Math::gaussian_quantile_c(0.5 * (1 - CL1), 1); // limit on deviation in sigmas

    cout << "Reject points with deviation more than " << setprecision(3) << devLimit << " sigma" << endl;

    // Determine lower point on momentum
    Double_t Plo = 0;
    while ((ref = (ParFileRef *)piter())) {
        Double_t P = ref->GetParameter(1);
        ref->GetObject("pangalp", pang);

        if (P < Plo)
            continue;

        if (pang->GetEntries() < 1000)
            Plo = P;
    }
    piter.Reset();
    cout << "Omit files with momentum lower than " << setprecision(5) << Plo << " GeV/c" << endl;

    TF1 *fang = new TF1("fang", "[0]+[1]*cos(0.01745329*x)+[2]*cos(2*0.01745329*x)", -180, 180);
    TF1 *fsig = new TF1("fsig",
                        "[0]+[1]*cos(0.01745329*x)+[2]*cos(2*0.01745329*x)+[3]*cos(3*0."
                        "01745329*x)+[4]*cos(4*0.01745329*x)+[5]*cos(5*0.01745329*x)",
                        -180, 180);

    const char *good = " (fit success)";
    const char *failed = " \033[41;1m(fit failed)\033[0m";

    Int_t rc1, rc2;

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 700);
    c1->Divide(1, 2);
    gStyle->SetOptStat(10);
    gStyle->SetOptFit(1011);

    ofstream fout("fitpar.dat");
    fout.precision(6);

    while ((ref = (ParFileRef *)piter())) {
        Double_t Theta = ref->GetParameter(0);
        Double_t P = ref->GetParameter(1);
        if (P < Plo)
            continue;

        cout.precision(4);
        cout << "Dip angle=" << Theta << " deg, momentum=" << P << " GeV/c" << endl;

        Double_t Beta = 1 / sqrt(1 + M * M / P / P);

        ref->GetObject("pangalp", pang);

        Double_t meanang = 0, meansig = 0;
        Int_t n = 0;
        for (Int_t bin = 1; bin <= nbins; bin++) {
            Double_t alpha = pang->GetBinCenter(bin);
            Double_t ang = pang->GetBinContent(bin);
            Double_t sig = pang->GetBinError(bin);
            Double_t ent = pang->GetBinEntries(bin);
            if (ent < 5) {
                cout << CT_FGREDBOLD << " rejected before fit: bin=" << bin << " alpha=" << alpha << " entries=" << ent
                     << CT_DEFAULT << endl;
                continue;
            }
            Double_t errang = sig / sqrt(ent);
            Double_t errsig = errang / sqrt(2);
            if (ang == 0 || sig == 0) {
                cout << CT_BGREDBOLD << " rejected before fit: bin=" << bin << " alpha=" << alpha << " angle=" << ang
                     << " sigma=" << sig << CT_DEFAULT << endl;
                gSystem->Sleep(1000);
                continue;
            }

            meanang += ang;
            meansig += sig;

            gang->SetPoint(n, alpha, ang);
            gang->SetPointError(n, 0, errang);

            gsig->SetPoint(n, alpha, sig);
            gsig->SetPointError(n, 0, errsig);

            n++;
        }

        meanang /= n;
        meansig /= n;

        if (n != gang->GetN()) {
            gang->Set(n);
            gsig->Set(n);
        }

        fang->SetParameters(0.5 * meanang, 0., 0.);
        if (Theta == 0) {
            fang->FixParameter(1, 0.);
            fang->FixParameter(2, 0.);
        } else {
            fang->ReleaseParameter(1);
            fang->ReleaseParameter(2);
        }

        fsig->SetParameters(0.5 * meansig, 0., 0., 0., 0., 0.);
        if (Theta == 0) {
            fsig->FixParameter(1, 0.);
            fsig->FixParameter(2, 0.);
            fsig->FixParameter(3, 0.);
            fsig->FixParameter(4, 0.);
            fsig->FixParameter(5, 0.);
        } else {
            fsig->ReleaseParameter(1);
            fsig->ReleaseParameter(2);
            fsig->ReleaseParameter(3);
            fsig->ReleaseParameter(4);
            fsig->ReleaseParameter(5);
        }

        gang->GetListOfFunctions()->Clear();
        gsig->GetListOfFunctions()->Clear();

        cout << " Fitting..." << endl;
        rc1 = gang->Fit(fang, "QR", "goff");
        rc2 = gsig->Fit(fsig, "QR", "goff");

        cout << "  Angle vs alpha fit (round 1): prob=" << (fang->GetProb() < 0.01 ? CT_FGREDBOLD : "")
             << fang->GetProb() << CT_DEFAULT << (rc1 ? failed : good) << "\n"
             << "  Error vs alpha fit (round 1): prob=" << (fsig->GetProb() < 0.01 ? CT_FGREDBOLD : "")
             << fsig->GetProb() << CT_DEFAULT << (rc2 ? failed : good) << endl;

        // Reject points with a large deviation from the fit
        for (Int_t i = 0; i < gang->GetN(); i++) {
            if (fabs(gang->GetY()[i] - fang->Eval(gang->GetX()[i])) > devLimit * gang->GetEY()[i]) {
                gang->RemovePoint(i);
                i--; // do not iterate point
            }
        }

        if (gang->GetN() < n) {
            cout << "  " << CT_FGREDBOLD << n - gang->GetN() << " angle points rejected" << CT_DEFAULT << "\n";
            if (n - gang->GetN() > 5)
                gSystem->Sleep(3000);
        }

        for (Int_t i = 0; i < gsig->GetN(); i++) {
            if (fabs(gsig->GetY()[i] - fsig->Eval(gsig->GetX()[i])) > devLimit * gsig->GetEY()[i]) {
                gsig->RemovePoint(i);
                i--; // do not iterate point
            }
        }

        if (gsig->GetN() < n) {
            cout << "  " << CT_FGREDBOLD << n - gsig->GetN() << " sigma points rejected" << CT_DEFAULT << "\n";
            if (n - gsig->GetN() > 5)
                gSystem->Sleep(3000);
        }

        // Fix badly defined parameters
        Bool_t fixAngPars = kFALSE;
        Bool_t fixSigPars = kFALSE;
        /*
                        for(Int_t ipar=1; ipar<fang->GetNpar(); ipar++) {
                                if( fabs(fang->GetParameter(ipar)) <
           2*fang->GetParError(ipar) ) { cout<<CT_FGCYANBOLD<<"  Angle fit
           par"<<ipar<<" big error: "
                                                <<100.*fang->GetParError(ipar)/fabs(fang->GetParameter(ipar))<<"%."
                                                <<CT_DEFAULT;
                                        fixAngPars=kTRUE;
                                }
                                if( fixAngPars ) {
                                        fang->FixParameter(ipar,0.);
                        cout<<CT_FGCYANBOLD<<"  Fix angle
           par"<<ipar<<CT_DEFAULT<<"\n";
                                }
                        }

                        for(Int_t ipar=1; ipar<fsig->GetNpar(); ipar++) {
                                if( fabs(fsig->GetParameter(ipar)) <
           2*fsig->GetParError(ipar) ) { cout<<CT_FGCYANBOLD<<"  Sigma fit
           par"<<ipar<<" big error: "
                                                <<100.*fsig->GetParError(ipar)/fabs(fsig->GetParameter(ipar))<<"%."
                                                <<CT_DEFAULT;
                                        fixSigPars=kTRUE;
                                }
                                if( fixSigPars ) {
                                        fsig->FixParameter(ipar,0.);
                        cout<<CT_FGCYANBOLD<<"  Fix sigma
           par"<<ipar<<CT_DEFAULT<<"\n";
                                }
                        }

                        for(Int_t ipar=fang->GetNpar()-1; ipar>0; ipar--) {
                                if( fabs(fang->GetParameter(ipar)) <
           2*fang->GetParError(ipar) ) { cout<<CT_FGCYANBOLD<<"  Angle fit
           par"<<ipar<<" big error: "
                                                <<100.*fang->GetParError(ipar)/fabs(fang->GetParameter(ipar))<<"%."
                                                <<CT_DEFAULT;
                        cout<<CT_FGCYANBOLD<<"  Fix angle
           par"<<ipar<<CT_DEFAULT<<"\n"; fixAngPars=kTRUE;
                                        fang->FixParameter(ipar,0.);
                                        rc1=gang->Fit(fang,"QR","goff");
                                } else
                        break;
                        }

                        for(Int_t ipar=fsig->GetNpar()-1; ipar>0; ipar--) {
                                if( fabs(fsig->GetParameter(ipar)) <
           2*fsig->GetParError(ipar) ) { cout<<CT_FGCYANBOLD<<"  Sigma fit
           par"<<ipar<<" big error: "
                                                <<100.*fsig->GetParError(ipar)/fabs(fsig->GetParameter(ipar))<<"%."
                                                <<CT_DEFAULT;
                        cout<<CT_FGCYANBOLD<<"  Fix sigma
           par"<<ipar<<CT_DEFAULT<<"\n"; fixSigPars=kTRUE;
                                        fsig->FixParameter(ipar,0.);
                                        rc2=gsig->Fit(fsig,"QR","goff");
                                } else
                                        break;
                        }
        */
        cout << " Refitting..." << endl;
        for (Int_t ipar = 1; ipar < fang->GetNpar(); ipar++) {
            fang->FixParameter(ipar, 0.);
        }
        for (Int_t ipar = 0; ipar < fang->GetNpar(); ipar++) {
            fang->ReleaseParameter(ipar);
            rc1 = gang->Fit(fang, "QR", "goff");
            if (fang->GetProb() > 0.1)
                break;
        }

        for (Int_t ipar = 1; ipar < fsig->GetNpar(); ipar++) {
            fsig->FixParameter(ipar, 0.);
        }
        for (Int_t ipar = 0; ipar < fsig->GetNpar(); ipar++) {
            fsig->ReleaseParameter(ipar);
            rc2 = gsig->Fit(fsig, "QR", "goff");
            if (fsig->GetProb() > 0.01)
                break;
        }

        /*		if( n>gang->GetN() || n>gsig->GetN() || fixAngPars || fixSigPars )
           { cout<<" Refitting..."<<endl; if( !fixAngPars )
                                        rc1=gang->Fit(fang,"QR","goff");
                                if( !fixSigPars )
                                        rc2=gsig->Fit(fsig,"QR","goff");
        */
        cout << "  Angle vs alpha fit (round 2): prob=" << (fang->GetProb() < 0.01 ? CT_FGREDBOLD : "")
             << fang->GetProb() << CT_DEFAULT << (rc1 ? failed : good) << "\n"
             << "  Error vs alpha fit (round 2): prob=" << (fsig->GetProb() < 0.01 ? CT_FGREDBOLD : "")
             << fsig->GetProb() << CT_DEFAULT << (rc2 ? failed : good) << endl;
        //		}

        if (rc1 || rc2)
            gSystem->Sleep(3000);

        cout.precision(4);
        cout << "  " << Beta << " " << Theta;
        fout << Beta << " " << Theta;
        cout << "  \033[1;33m";
        for (Int_t ipar = 0; ipar < fang->GetNpar(); ipar++) {
            cout << " " << fang->GetParameter(ipar);
            fout << " " << fang->GetParameter(ipar);
        }
        cout << "  \033[1;32m";
        for (Int_t ipar = 0; ipar < fsig->GetNpar(); ipar++) {
            cout << " " << fsig->GetParameter(ipar);
            fout << " " << fsig->GetParameter(ipar);
        }
        cout << CT_DEFAULT << endl;
        fout << "\n";

        c1->cd(1);
        gang->Draw("ape");
        gPad->Update();

        c1->cd(2);
        gsig->Draw("ape");
        gPad->Update();

        //		if( fang->GetProb()<0.01 || fsig->GetProb()<0.01 )
        // gSystem->Sleep(1000);
        if (fang->GetProb() < 0.01)
            gSystem->Sleep(1000);
    }

    fout.close();

    flist->Delete();
    delete flist;
}
