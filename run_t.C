{
	gROOT->Reset();
	gROOT->LoadMacro("getwidth.C");
	gROOT->LoadMacro("richres5.C");

        Bool_t perRing=kTRUE;
        const Float_t D=100.;
        const Float_t pixel=1.6;
        const Float_t eff=0.6;
	const Int_t nt=11;
	Double_t P=4000, n=1.07, t_tab[]={3.,5.,7.5,10.,12.5,15.,17.5,20.,25.,30.,35.};
	const char const *pmt = "burle";
	Double_t sRmr[3][nt], sumR[nt], sigMax=0;

	for(int i=0; i<nt; i++) {
		getres(t_tab[i],pixel,n,pmt,P,D-t_tab[i],eff);
                sRmr[0][i] = r.sRmr[0];
                if( perRing ) sRmr[0][i]/=sqrt(r.Npe[0]);
		sRmr[1][i] = r.sRmr[1];
                if( perRing ) sRmr[1][i]/=sqrt(r.Npe[0]);
		sRmr[2][i] = r.sRmr[2];
                if( perRing ) sRmr[2][i]/=sqrt(r.Npe[0]);
		sumR[i]=sqrt(sRmr[0][i]**2+sRmr[1][i]**2+sRmr[2][i]**2);
		if (sigMax<sumR[i]) sigMax=sumR[i];
	}

	TCanvas *c1 = new TCanvas("c1","Root canvas");
	c1->SetGrid();

	gStyle->SetLineWidth(2);
	TGraph gSig(nt,t_tab,sumR);
	gSig.SetLineStyle(1);
	gSig.SetLineWidth(3);
	TGraph gSigT(nt,t_tab,sRmr[0]);
	gSigT.SetMarkerStyle(1);
	gSigT.SetLineStyle(9);
	gSigT.SetLineColor(kBlue);
	gSigT.SetLineWidth(3);
	TGraph gSigP(nt,t_tab,sRmr[1]);
	gSigP.SetMarkerStyle(1);
	gSigP.SetLineStyle(2);
	gSigP.SetLineColor(kRed);
	gSigP.SetLineWidth(3);
	TGraph gSigD(nt,t_tab,sRmr[2]);
	gSigD.SetMarkerStyle(1);
	gSigD.SetLineColor(kMagenta);
	gSigD.SetLineStyle(5);
	gSigD.SetLineWidth(3);

        TString title=";thickness, mm;#sigma_{#theta} per ";
        if( perRing )
            title+="ring";
        else
            title+="photon";
        title+=", mrad";
        c1->DrawFrame(0,0,t_tab[nt-1]+5,TMath::Ceil(1.1*sigMax),title);
	gSig.Draw("c");
	gSigD.Draw("c");
	gSigT.Draw("c");
	gSigP.Draw("c");
}
