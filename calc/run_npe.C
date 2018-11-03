{
	gROOT->Reset();
	gROOT->LoadMacro("dataread.C");
	gROOT->LoadMacro("getwidth.C");
	gROOT->LoadMacro("richres4.C");
	const Int_t nt=100;
	Double_t maxt=300, t_tab[nt], n=1.07, p=4000;
	const char const *pmt = "hambs";
	Double_t npe[nt];

	for(int it=0; it<nt; it++) {
		t_tab[it]=maxt/nt*(it+1);
		getres(t_tab[it],1,n,pmt,p,100,.5,false);
		npe[it] = r.Npe[0];
	}

	TCanvas *c1 = new TCanvas("c1","Root canvas");
	c1->SetGrid();

	TGraph gnpe_t(nt,t_tab,npe);
	TH2F h2("h2","Number of photoelectrons vs thickness",
		2,0,maxt*1.05,2,0,npe[nt-1]*1.05);
	h2.SetXTitle("radiator thickness, mm");
	h2.SetYTitle("N_{pe}");
	h2.SetStats(0);
	h2.Draw();
	gnpe_t.Draw("c");
}
