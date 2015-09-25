#include "Spectra.h"


Spectra::Spectra()
{
    NP = 30;
    p1 = 0.04;
    p2 = 10.0;
    DM = 0;
    SM = 3;
    zeta = 0.05;

    if (DM==0) P = logspace(p1, p2, NP);
    else if (DM==1) P = linspace(p1, p2, NP);

    SPA = zeros(NP);
    SPV = zeros(NP);
    SPD = zeros(NP);
    SPE = zeros(NP);
    SPT = zeros(NP);
    SPI = zeros<int>(NP);
}

Spectra::Spectra(int np, double P1, double P2, int dm, int sm, bool pseudo, double Zeta)
{
    NP = np;
    p1 = P1;
    p2 = P2;
    DM = dm;
    SM = sm+pseudo*3;
    zeta = Zeta;

    if (DM==0) P = logspace(p1, p2, NP);
    else if (DM==1) P = linspace(p1, p2, NP);
    else if (DM==2) {
        if (p1>=1.0) P = linspace(p1, p2, NP);
        else {
            int Nshort = (int)(NP*0.5);
            int Nlong = NP-Nshort+1;
            double *Pshort = logspace(p1, 1.0, Nshort);
            double *Plong = linspace(1.0, p2, Nlong);
			P = new double[NP];
            for (int i = 0; i < Nshort; i++) P[i] = Pshort[i];
            for (int i = 1; i < Nlong; i++) P[i+Nshort-1] = Plong[i];
        }
    }

    SPA = zeros(NP);
    SPV = zeros(NP);
    SPD = zeros(NP);
    SPE = zeros(NP);
    SPT = zeros(NP);
    SPI = zeros<int>(NP);
}

Spectra::Spectra(int np, double *p, int sm, bool pseudo, double Zeta)
{
	NP = np;
    p1 = 0.04;
    p2 = 10.0;
    DM = 2;
	SM = sm + pseudo * 3;
	zeta = Zeta;

    P = p;

	SPA = zeros(NP);
	SPV = zeros(NP);
	SPD = zeros(NP);
	SPE = zeros(NP);
	SPT = zeros(NP);
	SPI = zeros<int>(NP);
}

void Spectra::reallocate()
{
    if ( 0 != P  )    delete [] P;
    if ( 0 != SPA)    delete [] SPA;
    if ( 0 != SPV)    delete [] SPV;
    if ( 0 != SPD)    delete [] SPD;
    if ( 0 != SPE)    delete [] SPE;
    if ( 0 != SPT)    delete [] SPT;
    if ( 0 != SPI)    delete [] SPI;

    SPA = zeros(NP);
    SPV = zeros(NP);
    SPD = zeros(NP);
    SPE = zeros(NP);
    SPT = zeros(NP);
    SPI = zeros<int>(NP);
}

double **Spectra::getSPData()
{
    double **data = new double*[6];
    for (int i = 0; i < 6; i++) data[i] = new double[NP];

    for (int i = 0; i < NP; i++)
    {
        data[0][i] =   P[i];
        data[1][i] = SPA[i];
        data[2][i] = SPV[i];
        data[3][i] = SPD[i];
        data[4][i] = SPE[i];
        data[5][i] = SPT[i];
    }
    return data;
}

void Spectra::setP(double *p, int np)
{
    DM = 2;

    if (NP!=np) {
        NP = np;
        if (0 != SPA )    delete [] SPA;
        if (0 != SPV )    delete [] SPV;
        if (0 != SPD )    delete [] SPD;
        if (0 != SPE )    delete [] SPE;
        if (0 != SPT )    delete [] SPT;
        if (0 != SPI )    delete [] SPI;

        SPA = zeros(NP);
        SPV = zeros(NP);
        SPD = zeros(NP);
        SPE = zeros(NP);
        SPT = zeros(NP);
        SPI = zeros<int>(NP);
    }

    if (0 != P) delete [] P;
    P = p;
}

void Spectra::setSPT(double *spt)
{
    for (int i = 0; i < NP; i++)
    {
        SPT[i] = spt[i];
    }
}

void Spectra::setSPT(double *p, double *spt, int np)
{
    DM = 2;

    if (NP!=np) {
        NP = np;
        if (0 != SPA )    delete [] SPA;
        if (0 != SPV )    delete [] SPV;
        if (0 != SPD )    delete [] SPD;
        if (0 != SPE )    delete [] SPE;
        if (0 != SPI )    delete [] SPI;

        SPA = zeros(NP);
        SPV = zeros(NP);
        SPD = zeros(NP);
        SPE = zeros(NP);
        SPT = zeros(NP);
        SPI = zeros<int>(NP);
    }

    if (0 !=   P) delete [] P;
    if (0 != SPT) delete [] SPT;
    P = p;
    SPT = spt;

    this->calc();
}

void Spectra::setSPT(double Tg, double PAF, double scale)
{
    double gamma = 0.9+(0.05-zeta)/(0.3+6.0*zeta);
    double eta1  = 0.02+(0.05-zeta)/(4.0+32.0*zeta);
    double eta2  = 1.0+(0.05-zeta)/(0.08+1.6*zeta);

    PAF = PAF*eta2;

    double Pi;
    for (int i = 0; i < NP; ++i)
    {
        Pi = P[i];
        if (Pi<0.1) SPT[i] = (1.0+(PAF-1.0)/0.1*Pi)*scale;
        else if (Pi <= Tg) SPT[i] = PAF*scale;
        else if (Pi <= 5.0*Tg) SPT[i] = PAF*pow((Tg/Pi),gamma)*scale;
        else SPT[i] = PAF*(pow(0.2,gamma)-eta1/eta2*(Pi-5.0*Tg))*scale;
    }
}

void Spectra::setAcc(double *Acc, int N, double DT)
{
    acc = Acc;
    n = N;
    dt = DT;
}

void Spectra::calc(bool all)
{
    if (all == false) spectrum(acc,&n,&dt,&zeta,P,&NP,SPA,SPI,&SM);
    else spectrumavd(acc,&n,&dt,&zeta,P,&NP,SPA,SPI,SPV,SPD,SPE,&SM);
}

double *Spectra::fitSP(double tol, int mit, int fm, double peak0)
{
    double *a = new double[n];
    peakScale(acc, n, peak0);
    fitspectrum(acc,&n,&dt,&zeta,P,&NP,SPT,a,&tol,&mit,&fm);
    return a;
}

void Spectra::fitError(double &Emax, double &Emean)
{
    Emax = 0.0;
    Emean = 0.0;

    double re;
    for (int i = 0; i < NP; i++) {
        re = fabs(SPA[i]-SPT[i])/SPT[i];
        if (re > Emax) Emax = re;
        Emean += re*re;
    }

    Emean = sqrt(Emean/NP);
}
