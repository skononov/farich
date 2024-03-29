#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
//#include <unistd.h>
#include <string>
#include <algorithm>
#include <set>

#include "TSystem.h"
#include "TTree.h"
#include "TRint.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"

#include "MLADescription.h"

using std::endl;
using std::cout;
using std::cerr;
using std::setprecision;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::vector;
using std::string;

static const char *optstring = "qn:s:N:D:T:p:b:B:e:o:O:i:m::a:C";
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
static bool polpar = false;
static int npol = 2;
static bool samethick = false;
static bool chromaticity = true;
static string absfn;
static string infn;
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
                << "/Rich/radiator/absDataFile" << (absfn.empty()?"none":absfn) << "\n"
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

static bool parse_param(string p)
{
    transform(p.begin(), p.end(), p.begin(), [](unsigned char c) { return tolower(c); });

    if (p.find_first_of("pol") == 0) {
        try {
            npol = stoi(p.substr(3));
        } catch (...) {
            return false;
        }
        if (npol < 1 || npol > 10)
            return false;
        polpar = true;
        if (p.back() == 's')
            samethick = true;
        else
            samethick = false;
        return true; 
    } else if (p == "nt") {
        polpar = false;
        return true;
    }

    return false;
}

static bool read_radiator_file(string fn,vector<pair<float,float>>& data)
{
    ifstream in(fn);
    if (in.fail()) return false;    
    
    float t, n;
    int line = 0; 
    while (1) {
        in >> t >> n;
        if (in.fail()) break;
        line++;
        if (t < 0) {
            cerr << fn << ": отрицательное значение толщины (" << t << ") в строке " << line << endl;
            return false;
        }
        if (n < 1.0) {
            cerr << fn << ": значение показателя (" << n << ") меньше 1.0 в строке " << line << endl;
            return false;
        }
        data.push_back({t,n});        
    }

    in.close();
        
    return true;
}

static void Usage(int status)
{
    cout
        << "Usage: " << progname << " [OPTIONS] qefile\n"
        << "Вычисляет разрешение детектора черенковских колец с многослойным фокусирующим аэрогелем.\n"
        << "Оптимизирует разрешение. Параметр qefile задает путь к файлу с данными о квантовой\n"
        << "эффективности фотонного детектора, представленых в двух колонках: длина волны (нм), эффективность (\%).\n\n"
        << " OPTIONS:\n"
        << "   -q                Задачный режим без графики и интерпретатора\n"
        << "   -e eff            Фактор к эффективности фотонного детектора, 0<eff<=1 (по умолчанию: " << efficiency << ")\n"
        << "   -p size           Размер квадратного пикселя фотонного детектора, мм (по умолчанию: " << pixelsize << ")\n"
        << "   -D distance       Расстояние от начала радиатора до фотонного детектора, мм (по умолчанию: " << D << ")\n"
        << "   -n ri1            Максимальный показатель преломления радиатора на 400 нм (по умолчанию: " << ri1 << ")\n"
        << "   -i filename       Задать описание радиатора из текстового файла с колонками:\n"
           "                       толщина слоя (мм), показатель преломления на 400 нм,\n"
           "                       начиная с первого слоя (по умолчанию: оптимизированный радиатор)\n"
        << "   -N nlayers        Число слоев аэрогеля, 0<nlayers<=100 (по умолчанию: " << nlayers << ")\n"
        << "   -T thickness      Толщина радиатора, мм (по умолчанию: " << T << ")\n"
        << "   -a filename       Задать спектр длины поглощения в аэрогеле из файла с колонками:\n"
           "                       длина волны (нм), длина поглощения (мм) (по умолчанию: нет поглощения)\n"
        << "   -s Lsc            Длина рассеяния на 400 нм, мм (по умолчанию: " << Lsc << ")\n"
        << "   -b beta           Оптимизировать радиатор для данной скорости частицы, 1/ri1<beta<=1 (по умолчанию: " << opbeta << ")\n"
        << "   -B beta           Сделать расчет для данной скорости частицы, 1/ri1<beta<=1 (по умолчанию:\n"
           "                       скорость частицы для оптимизации)\n"
        << "   -m[param]         Оптимизировать радиатор по угловой ошибке на трек, используя параметризацию param:\n"
           "                       nt, pol<k>[s] (k=1..10, s - одинаковая толщина слоев) (по умолчанию: nt)\n"
        << "   -C                Отключить дисперсию показателя преломления от длины волны (по умолчанию - включена)\n"
        << "   -o filename       Сохранить гистограмму распределения по радиусу в заданный ROOT-файл\n"
           "                       (по умолчанию: " << outfn << ")\n"
        << "   -O filename       Сохранить описание детектора в заданный командный файл Geant4.\n"
           "                       (по умолчанию: не cохранять).\n"
        << endl;
    exit(status);
}

#define PERC(a) 100. * (res1.a / res.a - 1)

int main(int argc, char *argv[])
{
    progname = argv[0];

    if (argc == 1)
        Usage(0);

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
        } else if (opt == 'i') {
            infn = optarg;
            if (gSystem->AccessPathName(infn.c_str())) { // файл не обнаружен
                cerr << infn << ": нет такого файла" << endl;
                return 1;
            }
        } else if (opt == 's') {
            float l = atof(optarg);
            if (l <= 0) {
                cerr << optarg << ": неправильно задана длина рассения" << endl;
                return 1;
            }
            Lsc = l;
        } else if (opt == 'a') {
            absfn = optarg;
            if (gSystem->AccessPathName(absfn.c_str())) { // файл не обнаружен
                cerr << absfn << ": нет такого файла" << endl;
                return 1;
            }
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
                if (!parse_param(optarg)) {
                    cerr << optarg << ": неизвестный метод параметризации" << endl;
                    return 1;
                }
            }
            minimize = true;
        } else if (opt == 'o') {
            outfn = optarg;
        } else if (opt == 'C') {
            chromaticity = false;
        } else if (opt == 'O') {
            macfn = optarg;
        } else {
            cerr << opt << ": неизвестная опция" << endl;
            Usage(1);
        }
    }
    //=========Проверка значений параметров программы==========//
    if (optind >= argc)
        Usage(1);
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
    if (minimize && nlayers == 1) {
        cerr << "Предупреждение: оптимизация однослойного радиатора невозможна" << endl;
        minimize = false;
    }
    
    qefn = argv[optind];
    if (gSystem->AccessPathName(qefn.c_str())) { // файл не обнаружен
        cerr << qefn << ": нет такого файла" << endl;
        return 1;
    }

    for(int iarg=0; iarg<argc; iarg++)
        cout << argv[iarg] << " ";
    cout << endl;

    //==========Расчет========//
    cout << "________________________________________________\n"
         << "Вычисление разрешения черенковского детектора.\n"
         << "Параметры:\n"
         << "  файл с квантовой эффективностью:              " << qefn << "\n"
         << "  фактор эффективности ФЭУ:                     " << efficiency << "\n"
         << "  размер пикселя ФЭУ:                           " << pixelsize << " мм\n"
         << "  расстояние от радиатора до ФЭУ:               " << D << " мм\n";
    if (infn.empty()) {
        cout << "  показатель преломления первого слоя на 400нм: " << ri1 << "\n"
             << "  число слоев радиатора:                        " << nlayers << "\n"
             << "  полная толщина радиатора:                     " << T << " мм\n";
    } else {
        cout << "  файл описания радиатора:                      " << infn << "\n";
    }
    cout << "  длина рассеяния в аэрогеле на 400 нм:         " << Lsc << " мм\n"
         << "  файл со спектром длины поглощения аэрогеле:   " << (absfn.empty()?"не задан":absfn) << "\n"
         << "  дисперсия показателя преломления аэрогеля:    " << (chromaticity?"включена":"отключена") << "\n"
         << "  оптимизация для скорости:                     " << opbeta << "\n"
         << "  расчет для скорости:                          " << beta << "\n";

    if (minimize) {
        if (polpar) {
            cout << "  Pol" << npol << "-оптимизация радиатора";
            if (samethick)
                cout << " с одинаковыми толщинами слоев\n";
            else
                cout << " с толщинами слоев из быстрой оптимизации\n";
        } else {
            cout << "  NT-оптимизация радиатора\n";
        }
    } else if (infn.empty()) {
        cout << "  только быстрая оптимизация радиатора\n";
    }
    if (!outfn.empty())
        cout << "  выходной root-файл: " << outfn << "\n";
    if (!macfn.empty())
        cout << "  файл макроса для Geant4: " << macfn << "\n";
    cout << "________________________________________________" << endl;

    TRint *app = 0;
    if (!batch) {
        int ac = 2;
        const char * const av[2] = {"rint", "-l"};
        app = new TRint("rint", &ac, const_cast<char**>(av));
    }

    //Чтение файла с данными квантовой эффективности
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
    cout.precision(6);

    //Чтение файла с данными длины поглощения света в аэрогеле
    Spectrum absLen;
    if (!absfn.empty()) {
        absLen.ReadFile(absfn.c_str());
        if (absLen.IsEmpty()) return 1;
        double awl1, awl2;
        absLen.GetRange(awl1,awl2);
        cout.precision(3);
        cout << "Длина поглощения света в аэрогеле определена от " << awl1 << " до " << awl2 << " нм, всего " << absLen.Size()
             << " значений." << endl;
        cout.precision(6);
        if (wl1 < awl1) { // add zero absorption length below the lowest wavelength
            absLen.AddEntry(awl1-10,0.);
            absLen.AddEntry(wl1,0.);
            cout << "Длина поглощения света продолжена в УФ до " << wl1 << " нм нулевыми значениями " << endl;
        }
    }

    double t0 = D-T; // proximity distance. may be not correct at this stage
    MLADescription mla(t0, opbeta);
    mla.SetScatteringLength(Lsc);
    mla.SetPDefficiency(phdeteff);
    mla.SetPixelSize(pixelsize);
    mla.SetAbsLength(absLen);
    mla.SetChromaticity(chromaticity);

    if (!infn.empty()) {
        vector<pair<float,float>> radiator;
        if (!read_radiator_file(infn, radiator))
            return 2;

        for(auto entry : radiator)
            mla.AddAlayer(entry.second,entry.first);

        T = mla.GetTotalThickness();
        t0 = D - T;
        if (t0 < 0.) {
            cerr << "Толщина радиатора T=" << T << " мм превышает расстояние между фотодетектором и аэрогелем D=" << D << " мм" << endl;
            return 2;
        }
        mla.SetProximityDistance(t0);
        nlayers = mla.GetNlayers();
        cout << "Заданный радиатор:" << endl;
    } else {
        double wl0 = mla.GetMaxSensitivityWL(T);
        cout << "Длина волны максимальной чувствительности к ЧИ: " << wl0 << " нм" << endl;
        mla.SetWavelength(wl0);

        mla.MakeFixed(nlayers, D, ri1);
        cout << "Радиатор после быстрой оптимизации:" << endl;
    }

    mla.Print("  ");

    struct MLADescription::Resolution res = mla.Calculate(beta, true, true);
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

    if (minimize) {
        bool opres;
        if (polpar)
            opres = mla.OptimizePol(nlayers, npol, D, ri1, samethick);
        else
            opres = mla.OptimizeNT(nlayers, D, ri1);

        if (opres) {
            cout << "Оптимизированный аэрогелевый радиатор" << endl;
            mla.Print("  ");

            struct MLADescription::Resolution res1 = mla.Calculate(beta, true, true);
            cout.precision(4);
            cout << "  Скорость частицы для расчета:                 " << res1.beta << "\n"
                 << "  Среднее число фотоэлектронов:                 " << res1.npe << " [" << setprecision(2) << PERC(npe) << "%]\n"
                 << "  Средний радиус:                               " << setprecision(4) << res1.radius << " мм [" << setprecision(2)
                 << PERC(radius) << "%]\n"
                 << "  Ошибка радиуса на 1 фотон (с учетом пикселя): " << setprecision(4) << res1.sigma1 << " мм ("
                 << res1.sigma1_px << " мм) [" << setprecision(2) << PERC(sigma1) << "% (" << PERC(sigma1_px) << "%)]\n"
                 << "  Ошибка радиуса на трек (с учетом пикселя):    " << setprecision(4) << res1.sigma_t << " мм ("
                 << res1.sigma_t_px << " мм) [" << setprecision(2) << PERC(sigma_t) << "% (" << PERC(sigma_t_px)
                 << "%)]\n"
                 << "  Ошибка угла на 1 фотон (с учетом пикселя):    " << setprecision(4) << res1.sigma1_ang
                 << " мрад (" << res1.sigma1_ang_px << " мрад) [" << setprecision(2) << PERC(sigma1_ang) << "% ("
                 << PERC(sigma1_ang_px) << "%)]\n"
                 << "  Ошибка угла на трек (с учетом пикселя):       " << setprecision(4) << res1.sigma_t_ang
                 << " мрад (" << res1.sigma_t_ang_px << " мрад) [" << setprecision(2) << PERC(sigma_t_ang) << "% ("
                 << PERC(sigma_t_ang_px) << "%)]" << endl;
            cout.precision(6);

            res = res1; // copy optimized radiator results
        } else
            cout << "Не удалось оптимизировать радиатор!" << endl;
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
