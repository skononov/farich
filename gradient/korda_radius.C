int Npoints=20;

void radius(Double_t nn, Double_t nk, Double_t La, Double_t Ltot, Double_t b)
{
    Double_t k = (nk-nn)/La;
    Double_t Lv = Ltot-La;
    Double_t b2 = b*b;

    if( TMath::Max(nn,nk) < 1/b ) {
        cerr<<"Particle is under Cherenkov threshold of radiator"<<endl;
        return;
    }

    TGraph* gR = new TGraph(Npoints);
    TString stan;
    stan.Form("1/sqrt((%6f+%6f*x)**2/([0]**2-1/%6f)-1)",nn,k,b2);
//    sRa.Form("sqrt((([0]*%4f)^2-1.)/((((%1f-%2f)/%3f*x+%2f)^2-[0]^2)*%4f^2+1.))",nk,nn,La,b);
//    x^2 => x**2 (Fortran notation)
    TF1* ftan = new TF1("ftan",stan,0.,La);

    for(int i=0;i<Npoints;i++)
    {
        Double_t x0 = La/(Npoints-1)*i;
        Double_t n0 = nn+k*x0;
        Double_t R = 0.;
        cout<<"R("<<x0<<") = ";
        if( n0 < 1/b ) { //below threshold at emission point
            cout<<"below threshold"<<endl;
            gR->SetPoint(i,x0,0.);
            continue;
        }
        if( n0*n0-1/b2 > 1 ) { //total internal reflection at radiator-air surface
            cout<<"total internal reflection"<<endl;
            gR->SetPoint(i,x0,0.);
            continue;
        }

        Double_t Rv = Lv/sqrt(1/(n0*n0-1/b2)-1);
        ftan->SetParameter(0,n0);
        Double_t Ra = ftan->Integral(x0,La);
        R = Ra + Rv;
        cout<<Ra<<" + "<<Rv<<" = "<<R<<endl;

        gR->SetPoint(i,x0,R);
    }

    gR->Draw("apc");
}
