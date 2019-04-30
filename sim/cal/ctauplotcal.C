#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"

Bool_t calLoaded = 0;
Float_t *Beta = 0, *Theta = 0, *a0 = 0, *a1 = 0, *a2 = 0, *b0 = 0, *b1 = 0, *b2 = 0, *b3 = 0, *b4 = 0, *b5 = 0;
int nBetaPoints = 0, nThetaPoints = 0;

void freecal()
{
    if (Beta)
        delete[] Beta;
    if (Theta)
        delete[] Theta;
    if (a0)
        delete[] a0;
    if (a1)
        delete[] a1;
    if (a2)
        delete[] a2;
    if (b0)
        delete[] b0;
    if (b1)
        delete[] b1;
    if (b2)
        delete[] b2;
    if (b3)
        delete[] b3;
    if (b4)
        delete[] b4;
    if (b5)
        delete[] b5;
    Beta = Theta = a0 = a1 = a2 = b0 = b1 = b2 = b3 = b4 = b5 = (Float_t *)0;
    nBetaPoints = nThetaPoints = 0;
    calLoaded = 0;
}

int loadcalfile(const char *fn)
{
    ifstream in(fn);
    if (in.fail()) {
        cerr << "!!! Can not open file " << fn << endl;
        return 0;
    }

    if (calLoaded) {
        cout << "Discard previously loaded calibration" << endl;
        freecal();
    }

    char str[256];
    int nlines = 0;

    while (1) {
        in.getline(str, 255);
        if (in.eof())
            break;
        nlines++;
    }
    in.clear();
    in.seekg(0, ios::beg);

    cout << "Calibration file " << fn << " has " << nlines << " lines" << endl;

    Beta = new Float_t[nlines];
    Theta = new Float_t[nlines];
    a0 = new Float_t[nlines];
    a1 = new Float_t[nlines];
    a2 = new Float_t[nlines];
    b0 = new Float_t[nlines];
    b1 = new Float_t[nlines];
    b2 = new Float_t[nlines];
    b3 = new Float_t[nlines];
    b4 = new Float_t[nlines];
    b5 = new Float_t[nlines];

    Bool_t error = 0;
    nBetaPoints = 0;
    nThetaPoints = 1;
    for (int i = 0; i < nlines; i++) {
        in >> Beta[i] >> Theta[i] >> a0[i] >> a1[i] >> a2[i] >> b0[i] >> b1[i] >> b2[i] >> b3[i] >> b4[i] >> b5[i];
        if (in.fail()) {
            cerr << "!!! Wrong format or early EOF. " << i << " lines read." << endl;
            error = 1;
            break;
        }

        if (i > 0 && Theta[i] != Theta[i - 1]) {
            if (Theta[i] < Theta[i - 1]) {
                cerr << "!!! Wrong order of theta values" << endl;
                error = 1;
            }
            nThetaPoints++;
        } else if (Theta[i] == Theta[0]) {
            if (i > 0 && Beta[i] < Beta[i - 1]) {
                cerr << "!!! Wrong order of beta values" << endl;
                error = 1;
            }
            nBetaPoints++;
        }
    }
    in.close();

    if (error) {
        freecal();
        return 0;
    }

    cout << nThetaPoints << " theta points, " << nBetaPoints << " beta points" << endl;

    if (nlines != nBetaPoints * nThetaPoints) {
        cerr << "!!! Inconsistency of beta and theta points with number of lines" << endl;
        freecal();
        return 0;
    }

    calLoaded = 1;
    return nlines;
}

void ctauplotcal(Float_t theta, const char *fn = "ctau_fitpar.dat")
{
    if (theta < 0 && theta > 40) {
        cerr << "!!! Illegal value of theta " << theta << endl;
        return;
    }

    if (!calLoaded)
        if (loadcalfile(fn) == 0)
            return;

    int ip = 0;
    Float_t *thetaPoints = new Float_t[nThetaPoints];
    for (int i = 0; i < nThetaPoints; i++) {
        thetaPoints[i] = Theta[i * nBetaPoints];
        if (thetaPoints[i] < theta)
            ip++;
    }

    if (ip == 0) {
        cout << "Theta value " << theta << " is lower than calibration limit " << thetaPoints[0] << endl;
        return;
    } else if (ip == nThetaPoints) {
        cout << "Theta value " << theta << " is higher than calibration limit " << thetaPoints[nThetaPoints - 1]
             << endl;
        return;
    }

    Double_t w1 = (thetaPoints[ip] - theta) / (thetaPoints[ip] - thetaPoints[ip - 1]),
             w2 = (theta - thetaPoints[ip - 1]) / (thetaPoints[ip] - thetaPoints[ip - 1]);
    cout << "Given theta is in the range (" << thetaPoints[ip - 1] << ", " << thetaPoints[ip] << ") degrees. "
         << "Weight1=" << w1 << " Weight2=" << w2 << endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 700);
    c1->DrawFrame(-200, 0, 200, 400,
                  TString::Format("Angle vs alpha calibration for dip "
                                  "angle=%.2f;#alpha, deg;#theta_{c}, mrad"));

    TCanvas *c2 = new TCanvas("c2", "c2", 500, 700);
    c2->DrawFrame(-200, 0, 200, 30,
                  TString::Format("Sigma vs alpha calibration for dip angle=%.2f;#alpha, "
                                  "deg;#sigma_{#theta_{c}}, mrad"));

    TF1 *fang = new TF1("fang", "[0]+[1]*cos(0.01745329*x)+[2]*cos(2*0.01745329*x)", -180, 180);
    TF1 *fsig = new TF1("fsig",
                        "[0]+[1]*cos(0.01745329*x)+[2]*cos(2*0.01745329*x)+[3]*cos(3*0."
                        "01745329*x)+[4]*cos(4*0.01745329*x)+[5]*cos(5*0.01745329*x)",
                        -180, 180);

    Float_t maxSig = 0;
    for (int i = 0; i < nBetaPoints; i++) {
        Int_t ip1 = (ip - 1) * nBetaPoints + i, ip2 = ip * nBetaPoints + i;

        fang->SetParameters(w1 * a0[ip1] + w2 * a0[ip2], w1 * a1[ip1] + w2 * a1[ip2], w1 * a2[ip1] + w2 * a2[ip2]);
        c1->cd();
        fang->DrawCopy("same");

        fsig->SetParameters(w1 * b0[ip1] + w2 * b0[ip2], w1 * b1[ip1] + w2 * b1[ip2], w1 * b2[ip1] + w2 * b2[ip2],
                            w1 * b3[ip1] + w2 * b3[ip2], w1 * b4[ip1] + w2 * b4[ip2], w1 * b5[ip1] + w2 * b5[ip2]);
        c2->cd();
        fsig->DrawCopy("same");

        if (maxSig < fsig->Eval(0))
            maxSig = fsig->Eval(0);
    }

    ((TH1F *)c2->GetPrimitive("hframe"))->SetMaximum(1.2 * maxSig);

    delete[] thetaPoints;
}
