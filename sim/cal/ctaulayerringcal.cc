namespace std {
}
using namespace std;

#include "ctaulayerringcal.C"
#include "TApplication.h"

int main(int argc, char *argv[])
{
    const char *fnpat = argc > 1 ? argv[1]
                                 : "/home/skononov/geant4/farich/runs/ctau/reco/"
                                   "pi_ctau_mppc2_px3_d200_mla4_*deg_*gev.root";
    const char *outfn = argc > 2 ? argv[2] : "fitpar_layers.dat";

    int ac = 3;
    char *av[] = {argv[0], "-l", "-x"};

    TApplication *app = new TApplication("TApplication", &ac, av);

    app->ExecuteFile("~/root/rootlogon.C");

    fitPDF = kTRUE;
    int rc = ctaulayerringcal(fnpat, outfn);

    delete app;

    return rc;
}
