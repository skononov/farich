#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <libgen.h>
#include <string>
#include <algorithm>
#include <set>

#include "TTree.h"
#include "TRint.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"

#include "MLADescription.h"

using namespace std;

static const char *optstring = "qn:s:N:D:T:p:b:B:e:o:O:m::";
static const char *progname;

//Параметры программы
static string qefn;
static bool batch = false;
static float ri1 = 1.07;
static float Lsc = 50.;
static int nlayers = 3;
static float D = 100.;
static float T = 25.;
static float efficiency = 1.0;
static float pixelsize = 0;
static float opbeta = 1.0;
static float beta = 0.;
static bool minimize = false;
static string param = "nt";
static bool polpar = false;
static int npol = 2;
static string outfn = "farichres.root";
static string macfn;

static void write_geant4_macfile(string macfn, MLADescription &mla)
{
    cout << "Сохраняем описание RICH в командный файл Geant4 " << macfn << endl;
    ofstream macfile(macfn.c_str());
    if (macfile.is_open()) {
        macfile << "/Rich/pmt/qeDataFile " << qefn << "\n"
                << "/Rich/pmt/detection " << efficiency << "\n"
                << "/Rich/radiator/scatterLength " << Lsc << " mm\n"
                << "/Rich/radiator/absDataFile none\n"
                << "/Rich/mode manual\n"
                << "/Rich/proximity " << mla.GetProximityDistance() << " mm\n";
        for (int l = 0; l < nlayers; l++) {
            macfile << "/Rich/radiator/addLayer " << setprecision(5) << mla.GetIndex(l) << " " << mla.GetThickness(l)
                    << " mm\n";
        }
        macfile << "/Rich/update" << endl;
    } else {
        cerr << macfn << ": не могу открыть файл на запись" << endl;
    }
}

string parse_param(string p, int &npol)
{
    transform(p.begin(), p.end(), p.begin(), [](unsigned char c) { return tolower(c); });
    if (p.find_first_of("pol") == 0) {
        try {
            npol = stoi(p.substr(3));
        } catch (...) {
            return string();
        }
        if (npol < 1 || npol > 10)
            return string();
        p = "pol";
    } else if (p != "nt")
        return string();

    return p;
}

void Usage(int status)
{
    cout
        << "Usage: " << progname << " [OPTIONS] qefile\n"
        << "Вычисляет разрешение детектора черенковских колец с многослойным фокусирующим аэрогелем.\n"
        << "Оптимизирует разрешение. Параметр qefile задает путь к файлу с данными о квантовой\n"
        << "эффективности фотонного детектора, представленых в двух колонках: длина волны (нм), эффективность (\%).\n\n"
        << " OPTIONS:\n"
        << "   -q                Задачный режим без графики и интерпретатора\n"
        << "   -e eff            Фактор к эффективности фотонного детектора, 0<eff<=1 (default: " << efficiency << ")\n"
        << "   -p size           Размер квадратного пикселя фотонного детектора, мм (default: " << pixelsize << ")\n"
        << "   -D distance       Расстояние от начала радиатора до фотонного детектора, мм (default: " << D << ")\n"
        << "   -n ri1            Максимальный показатель преломления радиатора на 400 нм (default: " << ri1 << ")\n"
        << "   -N nlayers        Число слоев аэрогеля, 0<nlayers<=100 (default: " << nlayers << ")\n"
        << "   -T thickness      Толщина радиатора, мм (default: " << T << ")\n"
        << "   -s Lsc            Длина рассеяния на 400 нм, мм (default: " << Lsc << ")\n"
        << "   -b beta           Оптимизировать радиатор для данной скорости частицы, 1/ri1<beta<=1 (default: "
        << opbeta << ")\n"
        << "   -B beta           Сделать расчет для данной скорости частицы, 1/ri1<beta<=1 (default: скорость частицы "
           "для оптимизации)\n"
        << "   -m [param]        Оптимизировать радиатор по угловой ошибке на трек, используя параметризацию param: "
           "nt, polN (N=1..10) (default: "
        << param << ")\n"
        << "   -o filename       Сохранить гистограмму распределения по радиусу в заданный ROOT-файл (default: "
        << outfn << ")\n"
        << "   -O filename       Сохранить описание детектора в заданный командный файл Geant4. По умолчанию - не "
           "сохранять.\n"
        << endl;
    exit(status);
}

#define PERC(a) 100. * (res1.a / res.a - 1)

int main(int argc, char *argv[])
{
    progname = argv[0];

    if (argc == 1)
        Usage(0);

    int npol;

    //=========Обработка параметров программы==========//
    int opt;
    while ((opt = getopt(argc, argv, optstring)) > 0) {
        if (opt == '?')
            Usage(1);
        else if (opt == 'q') {
            batch = true;
        } else if (opt == 'n') {
            float n = atof(optarg);
            if (n < 1 || n > 1.41) {
                cerr << optarg << ": сомнительный показатель преломления" << endl;
                return 1;
            }
            ri1 = n;
        } else if (opt == 's') {
            float l = atof(optarg);
            if (l <= 0) {
                cerr << optarg << ": неправильно задана длина рассения" << endl;
                return 1;
            }
            Lsc = l;
        } else if (opt == 'N') {
            int n = atoi(optarg);
            if (n < 1) {
                cerr << optarg << ": неправильно задано число слоев" << endl;
                return 1;
            }
            nlayers = n;
        } else if (opt == 'D') {
            float d = atof(optarg);
            if (d <= 0.0) {
                cerr << optarg << ": неправильно задано расстояние от радиатора до детектора" << endl;
                return 1;
            }
            D = d;
        } else if (opt == 'T') {
            float t = atof(optarg);
            if (t <= 0.0) {
                cerr << optarg << ": неправильно задана толщина радиатора" << endl;
                return 1;
            }
            T = t;
        } else if (opt == 'b') {
            float b = atof(optarg);
            if (b <= 0.0 || b > 1.0) {
                cerr << optarg << ": неправильно задана скорость частицы для оптимизации радиатора" << endl;
                return 1;
            }
            opbeta = b;
        } else if (opt == 'B') {
            float b = atof(optarg);
            if (b <= 0.0 || b > 1.0) {
                cerr << optarg << ": неправильно задана скорость частицы для расчета" << endl;
                return 1;
            }
            beta = b;
        } else if (opt == 'e') {
            float e = atof(optarg);
            if (e <= 0.0 || e > 1.0) {
                cerr << optarg << ": неправильно задан фактор эффективности фотонного детектора" << endl;
                return 1;
            }
            efficiency = e;
        } else if (opt == 'p') {
            float s = atof(optarg);
            if (s <= 0.0) {
                cerr << optarg << ": неправильно задан размер пикселя" << endl;
                return 1;
            }
            pixelsize = s;
        } else if (opt == 'm') {
            if (optarg) {
                param = parse_param(optarg, npol);
                if (param.empty()) {
                    cerr << optarg << ": неизвестный метод параметризации" << endl;
                    return 1;
                }
            }
            minimize = true;
        } else if (opt == 'o') {
            outfn = optarg;
        } else if (opt == 'O') {
            macfn = optarg;
        } else {
            cerr << opt << ": неизвестная опция" << endl;
            Usage(1);
        }
    }
    if (T >= D) {
        cerr << "T=" << T << ", D=" << D << ": толщина радиатора больше расстояния до фотонного детектора" << endl;
        return 1;
    }
    if (ri1 * opbeta <= 1.0) {
        cerr << "ri1=" << ri1 << " beta=" << opbeta << ": допороговая скорость для оптимизации в первом слое радиатора"
             << endl;
        return 1;
    }
    if (beta <= 0.)
        beta = opbeta;
    else if (ri1 * beta <= 1.0) {
        cerr << "ri1=" << ri1 << " beta=" << beta << ": допороговая скорость для расчета в первом слое радиатора"
             << endl;
        return 1;
    }

    if (optind >= argc)
        Usage(1);

    qefn = argv[optind];

    //==========Расчет========//
    cout << "________________________________________________\n"
         << "Вычисление разрешения черенковского детектора.\n"
         << "Параметры:\n"
         << "  файл с квантовой эффективностью:\n   " << qefn << "\n"
         << "  фактор эффективности ФЭУ:                     " << efficiency << "\n"
         << "  размер пикселя ФЭУ:                           " << pixelsize << " мм\n"
         << "  расстояние от радиатора до ФЭУ:               " << D << " мм\n"
         << "  показатель преломления первого слоя на 400нм: " << ri1 << "\n"
         << "  число слоев радиатора:                        " << nlayers << "\n"
         << "  полная толщина радиатора:                     " << T << " мм\n"
         << "  длина рассеяния в аэрогеле на 400 нм:         " << Lsc << " мм\n"
         << "  оптимизация для скорости:                     " << opbeta << "\n"
         << "  расчет для скорости:                          " << beta << "\n";
    if (minimize) {
        if (param == "nt")
            cout << "  оптимизация радиатора с методом параметризации NT\n";
        else
            cout << "  оптимизация радиатора с методом параметризации pol" << npol << "\n";
    }
    if (!outfn.empty())
        cout << "  выходной root-файл: " << outfn << "\n";
    if (!macfn.empty())
        cout << "  файл макроса для Geant4: " << macfn << "\n";
    cout << "________________________________________________" << endl;

    TRint *app = 0;
    if (!batch) {
        int ac = 2;
        char *av[2] = {"rint", "-l"};
        app = new TRint("rint", &ac, av);
    }

    Spectrum phdeteff(qefn.c_str());
    if (phdeteff.IsEmpty())
        return 1;

    if (efficiency < 1.0)
        phdeteff.Scale(efficiency);

    double wl1, wl2;
    phdeteff.GetRange(wl1, wl2);

    cout.precision(3);
    cout << "Квантовая эффективность определена от " << wl1 << " до " << wl2 << " нм, всего " << phdeteff.Size()
         << " значений." << endl;

    //Создаем многослойный аэрогелевый радиатор по идеальной модели без дисперсии для
    //длины волны максимальной чувствительности
    double t0 = D - T; // proximity distance
    MLADescription mla(t0, opbeta);
    mla.SetScatteringLength(Lsc);
    mla.SetPDefficiency(phdeteff);
    mla.SetPixelSize(pixelsize);

    mla.MakeFixed(nlayers, D, ri1);

    double wl0 = mla.GetMaxSensitivityWL(T);
    cout << "Длина волны максимальной чувствительности к ЧИ: " << wl0 << " нм" << endl;
    mla.SetWavelength(wl0);

    cout << "Исходный аэрогелевый радиатор:" << endl;
    mla.Print("  ");

    struct MLADescription::Resolution res = mla.Calculate(beta, true);
    cout.precision(4);
    cout << "  Скорость частицы для расчета:                 " << res.beta << "\n"
         << "  Среднее число фотоэлектронов:                 " << res.npe << "\n"
         << "  Средний радиус:                               " << res.radius << " мм\n"
         << "  Ошибка радиуса на 1 фотон (с учетом пикселя): " << res.sigma1 << " мм (" << res.sigma1_px << " мм)\n"
         << "  Ошибка радиуса на трек (с учетом пикселя):    " << res.sigma_t << " мм (" << res.sigma_t_px << " мм)\n"
         << "  Ошибка угла на 1 фотон (с учетом пикселя):    " << res.sigma1_ang << " мрад (" << res.sigma1_ang_px
         << " мрад)\n"
         << "  Ошибка угла на трек (с учетом пикселя):       " << res.sigma_t_ang << " мрад (" << res.sigma_t_ang_px
         << " мрад)" << endl;
    cout.precision(6);

    if (minimize && nlayers > 1) {
        bool opres;
        if (param == "pol")
            opres = mla.OptimizePol(nlayers, npol, D, ri1);
        else
            opres = mla.OptimizeNT(nlayers, D, ri1);

        if (opres) {
            cout << "Оптимизированный аэрогелевый радиатор" << endl;
            mla.Print("  ");

            struct MLADescription::Resolution res1 = mla.Calculate(beta, true);
            cout.precision(4);
            cout << "  Скорость частицы для расчета:                 " << res1.beta << "\n"
                 << "  Среднее число фотоэлектронов:                 " << res1.npe << "\n"
                 << "  Средний радиус:                               " << res1.radius << " мм [" << setprecision(2)
                 << PERC(radius) << "%]\n"
                 << "  Ошибка радиуса на 1 фотон (с учетом пикселя): " << setprecision(4) << res1.sigma1 << " мм ("
                 << res.sigma1_px << " мм) [" << setprecision(2) << PERC(sigma1) << "% (" << PERC(sigma1_px) << "%)]\n"
                 << "  Ошибка радиуса на трек (с учетом пикселя):    " << setprecision(4) << res1.sigma_t << " мм ("
                 << res.sigma_t_px << " мм) [" << setprecision(2) << PERC(sigma_t) << "% (" << PERC(sigma_t_px)
                 << "%)]\n"
                 << "  Ошибка угла на 1 фотон (с учетом пикселя):    " << setprecision(4) << res1.sigma1_ang
                 << " мрад (" << res.sigma1_ang_px << " мрад) [" << setprecision(2) << PERC(sigma1_ang) << "% ("
                 << PERC(sigma1_ang_px) << "%)]\n"
                 << "  Ошибка угла на трек (с учетом пикселя):       " << setprecision(4) << res1.sigma_t_ang
                 << " мрад (" << res.sigma_t_ang_px << " мрад) [" << setprecision(2) << PERC(sigma_t_ang) << "% ("
                 << PERC(sigma_t_ang_px) << "%)]" << endl;
            cout.precision(6);

            res = res1; // copy optimized radiator results
        } else
            cout << "Не удалось оптимизировать радиатор!" << endl;
    } else if (minimize && nlayers == 1) {
        cerr << "Оптимизация однослойного радиатора невозможна!" << endl;
    }

    TFile *outfile = nullptr;
    if (!outfn.empty()) {
        cout << "Сохраняем распределение фотонов по радиусу в " << outfn << endl;
        outfile = new TFile(outfn.c_str(), "RECREATE");
        if (!outfile->IsOpen()) {
            cerr << outfn << ": не могу открыть файл на запись" << endl;
            delete outfile;
            outfile = nullptr;
        }
    }

    TH1D hrad("hrad", "Radius photoelectron distribution;radius, mm;dNpe/dR, mm^{-1}", MLADescription::Nr, res.rmin,
              res.rmax);
    for (size_t i = 0; i < res.s.size(); i++) {
        hrad.SetBinContent(i + 1, res.s[i]);
    }

    double *xbins = new double[nlayers + 1];
    xbins[0] = 0.;
    for (int i = 0; i < nlayers; i++)
        xbins[i + 1] = xbins[i] + mla.GetThickness(i);
    TH1D hri("hri", "Refractive index profile;radius, mm;refractive index", nlayers, xbins);
    for (int bin = 1; bin <= nlayers; bin++)
        hri.SetBinContent(bin, mla.GetIndex(bin - 1));

    if (outfile) {
        MLADescription::CalcTuple tuple;
        TTree *T = new TTree("T", "FARICH resolution detailed calculation results");
        T->Branch("data", &tuple, "l/I:r/F:wl/F:x0/F:s/F");
        T->SetMarkerStyle(7);

        for (MLADescription::CalcTuple &_tuple : res.data) {
            tuple = _tuple;
            T->Fill();
        }

        // Writing ROOT file
        hrad.Write();
        hri.Write();
        T->Write();
    }

    if (!macfn.empty())
        write_geant4_macfile(macfn, mla);

    if (!batch) {
        if (outfile)
            outfile->ReOpen("READ");

        TCanvas *c1 = new TCanvas("c1", "c1", 802, 428);
        c1->Divide(2, 1);
        c1->cd(1);
        hrad.Draw("l");
        c1->cd(2);
        hri.Draw();

        app->Run(kTRUE);
    }

    if (outfile)
        delete outfile;

    return 0;
}
