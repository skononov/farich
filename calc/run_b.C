{
	gROOT->Reset();
	gROOT->LoadMacro("dataread.C");
	gROOT->LoadMacro("getwidth.C");
	gROOT->LoadMacro("richres4.C");
	const Int_t np=10;
	Double_t Pmin=1.35, Pmax=3, p_tab[np], n=1.07;
	const char const *pmt = "gaas_intevac";
	Double_t sR[2][np], sigMax=0;


	for(int ip=0; ip<np; ip++) {
		p_tab[ip]=Pmin + (Pmax-Pmin)/(np-1)*ip;
		getres(12,1,n,pmt,1e3*p_tab[ip],88,.5,false);
		sR[0][ip] = r.sR[1][0];
		sR[1][ip] = r.sR[1][2];
		if (sigMax<TMath::Max(r.sR[1][0],r.sR[1][2]))
			sigMax=TMath::Max(r.sR[1][0],r.sR[1][2]);
	}

	TCanvas *c1 = new TCanvas("c1","Root canvas");
	c1->SetGrid();

	TGraph gSigT(np,p_tab,sR[0]);
	TGraph gSigD(np,p_tab,sR[1]);
	gSigT.SetMarkerStyle(20);
	gSigD.SetMarkerStyle(24);

	TH2F h2("h2","Radius resolutions vs momentum",
		2,Pmin*.98,Pmax*1.02,2,0,1.1*sigMax);
	h2.SetXTitle("momentum, GeV/c");
	h2.SetYTitle("#sigma_{r}, mm");
	h2.SetStats(0);
	h2.Draw();
	gSigT.Draw("lp");
	gSigD.Draw("lp");
}
