
Double_t tan_vs_x(Double_t *xx, Double_t *par)
{
    Double_t x = *xx;
    Double_t x0 = par[0], n1 = par[1], n2 = par[2], nt = par[3];
    Double_t n0 = n1*(1-x0) + n2*x0;
    if( nt>n0 ) return 0.;

    Double_t n;
    if( x<0. || x>1. )
        n = 1.;
    else
        n = n1*(1-x) + n2*x;

    return 1./sqrt(n*n/(n0*n0-nt*nt)-1);
}

TF1* ftanx = new TF1("ftanx",tan_vs_x,0.,1.,4);

Double_t radius_vs_origin_int(Double_t *xx, Double_t* par)
{
    Double_t x0 = *xx;
    //                      n1       n2      nt
    Double_t pp[4] = { x0, par[0], par[1], par[2] };
    Double_t Lv = par[3]; //proximity gap divided by radiator thickness

    ftanx->SetParameters(pp);

    return ftanx->Eval(1.1)*Lv + ftanx->Integral(x0,1.);
}

Double_t radius_vs_origin_ana(Double_t *xx, Double_t* par)
{
    Double_t x0 = *xx;
    Double_t n1 = par[0], n2 = par[1], nt = par[2];
    Double_t Lv = par[3]; //proximity gap divided by radiator thickness
    Double_t n0 = n1*(1-x0)+n2*x0;
    if( n0<nt ) return 0.;

    Double_t nn0t = n0*n0-nt*nt;

    return sqrt(nn0t)/(n2-n1)*log((n2+sqrt(n2*n2-nn0t))/(n0+nt)) + Lv/sqrt(1/nn0t-1);
}

TF1* frado = new TF1("frado",radius_vs_origin_ana,0.,1.,4);

void plot_radius_vs_origin(Double_t T, Double_t D, Double_t n1, Double_t n2, Double_t beta)
{
    if( D<T ) {
        cerr<<"Dimension less than thickness!"<<endl;
        return;
    }
    Double_t Lv = D/T-1;

    if( n1<1. || n2<1. ) {
        cerr<<"Illegal refractive index values!"<<endl;
        return;
    }
    Double_t nt = 1/beta;

    if( TMath::Max(n1,n2)<nt ) {
        cerr<<"Refractive index everywhere is less than threshold value "<<nt<<endl;
        return;
    }

    frado->SetParameters(n1,n2,nt,Lv);

    frado->Draw();
}

void lineargrad(Double_t T, Double_t D, Double_t n2, Double_t beta)
{
    if( D<T ) {
        cerr<<"Dimension less than thickness!"<<endl;
        return;
    }
    Double_t Lv = D/T-1;

    if( n2<1. ) {
        cerr<<"Illegal refractive index value!"<<endl;
        return;
    }
    Double_t nt = 1/beta;

    if( n2<nt ) {
        cerr<<"Refractive index is less than threshold value "<<nt<<endl;
        return;
    }

    //initial estimate for n1
    Double_t n1 = nt+(n2-nt)*TMath::Power(Lv/(1+Lv),2);

    frado->SetParameters(n1,n2,nt,Lv);

    //goal radius
    Double_t Rg = frado->Eval(1.);

    Double_t R0 = frado->Eval(0.);

    while( fabs(R0/Rg-1)>1e-12 ) {
        Double_t A = R0/sqrt(n1*n1-nt*nt);
        n1 = sqrt(nt*nt+Rg*Rg/(A*A));
        frado->SetParameter(0,n1);
        R0 = frado->Eval(0.);
    }

    cout<<"Found n1 = "<<n1<<endl;

    plot_radius_vs_origin(T,D,n1,n2,beta);
}
