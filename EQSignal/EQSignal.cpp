#include "EQSignal.h"
#include "eqs.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#define NS 2

using namespace std;

static inline int countFileLines(const char *filename)
{
	int n = 0;
	ifstream in(filename, ios::in);
	string line;

	while (getline(in, line))
	{
		++n;
	}

	in.close();

	return n;

}

EQSignal::EQSignal()
{
    n = 4096;
    nsp = 2;
    dt = 0.02;
    v0 = 0.0;
    d0 = 0.0;

	t = linspace(0.0, dt*(n - 1), n);
    acc0 = zeros(n);
    acc  = zeros(n);
    vel  = zeros(n);
    disp = zeros(n);
    td   = zeros(n);
    tv   = zeros(n);
    ta   = zeros(n);
    rd   = zeros(n);
    rv   = zeros(n);
    ra   = zeros(n);
    Ia   = zeros(n);
    Iv   = zeros(n);
    Id   = zeros(n);

	res.ra = ra;
	res.rv = rv;
	res.rd = rd;
	res.ku = zeros(n);

    nfft = nextpow2(n);
    npsd = nfft/NS;
    af = new std::complex<double>[nfft];
    double fn = 0.5/dt;
    freqs = linspace(0.0,fn,nfft/2+1);
    ampf = zeros(nfft/2+1);
    angf = zeros(nfft/2+1);
    dangf = zeros(nfft/2+1);
    fpsd = linspace(0.0,fn,npsd/2+1);
    psd = zeros(npsd/2+1);

	sp = new Spectra[nsp];
    zeta = new double[nsp];
}

EQSignal::EQSignal(int N, double DT, double V0, double D0)
{
    n = N;
    nsp = 2;
	dt = DT;
    v0 = V0;
    d0 = D0;

	t = linspace(0.0, dt*(n - 1), n);
    acc0 = zeros(n);
    acc  = zeros(n);
    vel  = zeros(n);
    disp = zeros(n);
    td   = zeros(n);
    tv   = zeros(n);
    ta   = zeros(n);
    rd   = zeros(n);
    rv   = zeros(n);
    ra   = zeros(n);
    Ia   = zeros(n);
    Iv   = zeros(n);
    Id   = zeros(n);

	res.ra = ra;
	res.rv = rv;
	res.rd = rd;
	res.ku = zeros(n);

    nfft = nextpow2(n);
    npsd = nfft/NS;
    af = new std::complex<double>[nfft];
    double fn = 0.5/dt;
    freqs = linspace(0.0,fn,nfft/2+1);
    ampf = zeros(nfft/2+1);
    angf = zeros(nfft/2+1);
    dangf = zeros(nfft/2+1);
    fpsd = linspace(0.0,fn,npsd/2+1);
    psd = zeros(npsd/2+1);

	sp = new Spectra[nsp];
    zeta = new double[nsp];
}

EQSignal::EQSignal(double *a, int N, double DT, double V0, double D0)
{
    n   = N;
    nsp = 2;
    dt  = DT;
    v0  = V0;
    d0  = D0;

    t    = linspace(0.0, dt*(n - 1), n);
    acc0 = zeros(n);
    acc  = zeros(n);
    vel  = zeros(n);
    disp = zeros(n);
    td   = zeros(n);
    tv   = zeros(n);
    ta   = zeros(n);
    rd   = zeros(n);
    rv   = zeros(n);
    ra   = zeros(n);
    Ia   = zeros(n);
    Iv   = zeros(n);
    Id   = zeros(n);

	res.ra = ra;
	res.rv = rv;
	res.rd = rd;
	res.ku = zeros(n);

    for (int i = 0; i < n; ++i)
    {
        t[i]    = i*dt;
        acc[i]  = a[i];
        acc0[i] = acc[i];
    }

    this->a2vd();

    nfft = nextpow2(n);
    npsd = nfft/NS;
    af = new std::complex<double>[nfft];
    double fn = 0.5/dt;
    freqs = linspace(0.0,fn,nfft/2+1);
    ampf = zeros(nfft/2+1);
    angf = zeros(nfft/2+1);
    dangf = zeros(nfft/2+1);
    fpsd = linspace(0.0,fn,npsd/2+1);
    psd = zeros(npsd/2+1);

	sp = new Spectra[nsp];
    zeta = new double[nsp];
}

EQSignal::~EQSignal()
{
    //if (t != 0)    delete [] t;
    //if (acc0 != 0) delete [] acc0;
    //if (acc != 0)  delete [] acc;
    //if (vel != 0)  delete [] vel;
    //if (disp != 0) delete [] disp;
    //if (Ia != 0)   delete [] Ia;
    //if (Iv != 0)   delete [] Iv;
    //if (Id != 0)   delete [] Id;
    //if (ta != 0)   delete [] ta;
    //if (tv != 0)   delete [] tv;
    //if (td != 0)   delete [] td;
    //if (ra != 0)   delete [] ra;
    //if (rv != 0)   delete [] rv;
    //if (rd != 0)   delete [] rd;
    //if (sp != 0)   delete[] sp;
    //if (zeta != 0) delete[] zeta;
    //if (af != 0)   delete [] af;
    //if (freqs != 0)   delete [] freqs;
    //if (ampf != 0)   delete [] ampf;
    //if (angf != 0)   delete [] angf;
    //if (fpsd != 0)   delete [] fpsd;
    //if (psd != 0)   delete [] psd;
}

void EQSignal::readtxt(const char *filename, double DT, bool NORM)
{
    n = countFileLines(filename);
    dt = DT;

    this->reallocateTH();

    ifstream in(filename);

    for (int i = 0; i < n; ++i)
    {
        t[i] = i*dt;
        in >> acc[i];
        acc0[i] = acc[i];
    }

    in.close();

    if (NORM) this->norm();

    this->a2vd();

}

void EQSignal::readtxt(QString filename, double DT, bool NORM)
{
    this->readtxt(filename.toUtf8().data(),DT,NORM);
}

void EQSignal::readnga(const char *filename, bool NORM)
{
    ifstream in(filename, ios::in);
    string line;
    for (int i = 0; i < 3; i++) getline(in, line);
	in >> n >> dt >> line >> line;

    this->reallocateTH();

    for (int i = 0; i < n; i++)
    {
        t[i] = i*dt;
        in >> acc[i];
		acc0[i] = acc[i];
    }

    in.close();

    if (NORM) this->norm();

    this->a2vd();
}

void EQSignal::readnga(QString filename, bool NORM)
{
    this->readnga(filename.toUtf8().data(),NORM);
}

void EQSignal::resample(int r)
{
    int nr = int((double)(n-1)/(double)r) + 1;
    double *ar = new double[nr];
    fftresample(acc,&n,&r,ar,&nr);

    n = nr;
    dt = dt*r;
    this->reallocateTH();

    for (int i = 0; i < n; ++i)
    {
        acc[i] = ar[i];
        acc0[i] = acc[i];
    }

    this->a2vd();

    delete [] ar;
}

void EQSignal::interpolate()
{

}

void EQSignal::recover()
{
    for (int i = 0; i < n; i++) acc[i] = acc0[i];
}

void EQSignal::trim(int ind1, int ind2)
{
    n = ind2 - ind1 + 1;
	double *a = arraySlice(acc, ind1, ind2);
    this->reallocateTH();

	for (int i = 0; i < n; i++) { acc[i] = a[i]; acc0[i] = a[i]; }
    this->a2vd();
}

void EQSignal::confirm()
{
    for (int i = 0; i < n; i++) acc0[i] = acc[i];
}

void EQSignal::a2vd(bool raw, bool rat)
{
    if (raw) this->recover();

    if (rat) ratacc2vd(acc,vel,disp,&n,&dt,&v0,&d0);
    else acc2vd(acc,vel,disp,&n,&dt,&v0,&d0);
}

void EQSignal::norm()
{
	normalize(acc, n);
	normalize(acc0, n);

    this->a2vd();
}

void EQSignal::calcAriasIntensity()
{
    ariasintensity(acc,Ia,&n);
    ariasintensity(vel,Iv,&n);
    ariasintensity(disp,Id,&n);
}

void EQSignal::copyAccFrom(EQSignal *eqs)
{
    double *a = eqs->getAcc();
    int N = eqs->getN();

    for (int i = 0; i < min(n,N); ++i) acc[i] = a[i];

    if (n>N) for (int i = N; i < n; ++i) acc[i] = 0.0;

}

int *EQSignal::autoTrimEdges(int method, double thd1, double thd2, bool EZ)
{
    int *edges = new int[2];
	edges[0] = 0;
	edges[1] = n - 1;
    double pk, th1, th2;
    switch (method) {
    case 0:
        this->calcAriasIntensity();
        for (int i = 0; i < n; ++i)
        {
            if (Ia[i]>thd1) { edges[0]=i; break; }
        }
        for (int i = 0; i < n; ++i)
        {
			if (Ia[n - 1 - i]<thd2) { edges[1] = n - 1 - i; break; }
        }
        break;
	case 1:
        pk = fabs(peak(acc, n));
        th1 = pk*thd1;
        th2 = pk*thd2;
		for (int i = 0; i < n; ++i)
		{
            if (fabs(acc[i])>th1) { edges[0] = i; break; }
		}
		for (int i = 0; i < n; ++i)
		{
            if (fabs(acc[n-1-i])<(pk-th2)) { edges[1] = n - 1 - i; break; }
		}
		break;
    default:
        break;
    }

    if (EZ && edges[0]>0) {
        while (acc[edges[0]]*acc[edges[0]-1]>0) edges[0]--;
    }
	if (EZ && edges[1]<n) {
		while (acc[edges[1]] * acc[edges[1] + 1]>0) edges[1]++;
	}
    return edges;
}

void EQSignal::detrend(int method, int oh, int ol, bool raw)
{
    if (raw) this->recover();

    if (method == 0) polyblc(acc,&n,&oh,&ol,&dt,&v0,&d0);
    else polydetrend(acc,&n,&oh,&ol);

    this->a2vd();
}

void EQSignal::align(int method, int ntp, int oh, int ol,
                     bool raw, bool EWZ)
{
    if (raw)
    {
        this->recover();
        this->a2vd();
    }

    for (int i = 0; i < n; i++)
    {
        td[i] = disp[i];
        tv[i] = vel[i];
        ta[i] = acc[i];
    }

    if (oh==0) {
        bilinearDetrend(td,n);
        bilinearDetrend(tv,n);
        bilinearDetrend(ta,n);
    }
    else
    {
        polydetrend(td,&n,&oh,&ol);
        polydetrend(tv,&n,&oh,&ol);
        polydetrend(ta,&n,&oh,&ol);
    }

    int *tp = new int[ntp];
    int pl = 2;
    int ph = pl + ntp - 1;

    if (ntp == 1) tp[0] = int(0.618*n);
    else for (int i = 0; i < ntp; i++) tp[i] = int((i+1)*1.0/ntp*(n));

    switch (method) {
    case 0:
		if (EWZ) { td[n - 1] = 0.0; }
        targetdc(acc,td,&n,tp,&ntp,&ph,&pl,&dt,&v0,&d0);
        break;
    case 1:
		if (EWZ) { td[n - 1] = 0.0; tv[n - 1] = 0.0; }
        targetdvc(acc,td,tv,&n,tp,&ntp,&ph,&pl,&dt,&v0,&d0);
        break;
    case 2:
		if (EWZ) { td[n - 1] = 0.0; tv[n - 1] = 0.0; ta[n - 1] = 0.0; }
        targetdvac(acc,td,tv,ta,&n,tp,&ntp,&ph,&pl,&dt,&v0,&d0);
        break;
    default:
        break;
    }

    delete [] tp;
    this->a2vd();
}

void EQSignal::endAlign(int ntp, bool raw, int IZC, bool AccOnly)
{
    if (raw)
    {
        recover();
        a2vd();
    }

    beginLinearDetrend(acc,n,IZC);
    endLinearDetrend(acc,n,IZC);
    a2vd();

    if (AccOnly) return;

    for (int i = 0; i < n; i++)
    {
        td[i] = disp[i];
        tv[i] = vel[i];
        ta[i] = acc[i];
    }

    endLinearDetrend(td,n,IZC);
    endLinearDetrend(tv,n,IZC);
//    endLinearDetrend(ta,n,IZC);

    int *tp = new int[ntp];
    int pl = 2;
    int ph = pl + ntp - 1;

    if (ntp == 1) tp[0] = n;
    else for (int i = 0; i < ntp; i++) tp[i] = int((i+1)*1.0/ntp*n);

    targetdvac(acc,td,tv,ta,&n,tp,&ntp,&ph,&pl,&dt,&v0,&d0);

    delete [] tp;
//    confirm();
    a2vd();
}

void EQSignal::filt(int ftype, int order, double f1, double f2, bool raw)
{

    if (raw) this->recover();

	bwfilt(acc, n, dt, ftype, order, f1, f2);

    this->a2vd();
}

void EQSignal::setupSP(int NSP, int np, double P1, double P2, int dm, int sm, bool pseudo, double *Zetas)
{
    if (sp != 0) delete [] sp;
    if (zeta != 0) delete [] zeta;

	nsp = NSP;
    sp = new Spectra[nsp];
	
    zeta = Zetas;

    for (int i = 0; i < NSP; i++)
        sp[i] = Spectra(np,P1,P2,dm,sm,pseudo,zeta[i]);
}

void EQSignal::setupSP(int NSP, int np, double *p, int sm, bool pseudo, double *Zetas)
{
    if (sp != 0) delete [] sp;
    if (zeta != 0) delete [] zeta;

    nsp = NSP;
    sp = new Spectra[nsp];

    zeta = Zetas;

    for (int i = 0; i < NSP; i++)
        sp[i] = Spectra(np,p,sm,pseudo,zeta[i]);
}

void EQSignal::calcSP(bool allSP)
{
    for (int i = 0; i < nsp; i++)
    {
        Spectra *spi = sp + i;
        spi->setAcc(acc,n,dt);
        spi->calc(allSP);
    }
}

void EQSignal::calcSP(int i, bool allSP)
{
    Spectra *spi = sp + i;
    spi->setAcc(acc,n,dt);
    spi->calc(allSP);
}

void EQSignal::calcFFT()
{
	double *a = zeros(nfft);
	for (int i = 0; i < n; ++i) a[i] = acc[i];

    fft(a,af,&nfft);

    for (int i=0; i<nfft/2+1; ++i) {
        ampf[i] = 2.0*abs(af[i])/((double)nfft);
        angf[i] = arg(af[i]);
        if (angf[i]<0.0) angf[i] = angf[i]+PI*2.0;
        if (i>0) {
            dangf[i] = angf[i] - angf[i-1];
            if (dangf[i]>0.0) dangf[i] = dangf[i]-PI*2.0;
            else if (dangf[i]<-PI*2.0) dangf[i] = dangf[i]+PI*2.0;
        }
    }
    delete [] a;
}

void EQSignal::calcPSD(double olr,bool win)
{
    double *p = new double[npsd];
    int w;
    if (win) w = 1;
    else w = -1;
    welch(acc,&n,&npsd,&olr,p,&w);
    for (int i=0; i<npsd/2+1; ++i) {
        psd[i] = 2.0*p[i];
    }
}

void EQSignal::setSPT(double Tg, double amax, double scale)
{
    for (int i = 0; i < nsp; i++) (sp+i)->setSPT(Tg,amax,scale);
}

void EQSignal::fitSP(int i, double tol, int mit, int fm, double peak0)
{
    acc = sp[i].fitSP(tol,mit,fm,peak0);
    this->confirm();
}



double EQSignal::getDR()
{
    double *d = new double[n];
    for (int i = 0; i < n; i++) d[i] = disp[i];

    int oh = 1;
    int ol = 0;
    polydetrend(d,&n,&oh,&ol);

    double md = mean(disp,n);
    double rd = rms(d,n);

    delete [] d;

    return fabs(md / rd);
}

double EQSignal::getAR()
{
    double *d = new double[n];
    for (int i = 0; i < n; i++) d[i] = disp[i];

    int oh = 2;
    int ol = 0;
    polydetrend(d,&n,&oh,&ol);

    double md = rms(disp,n);
    double rd = rms(d,n);

    delete [] d;

    return fabs(md / rd)-1.0;
}

double *EQSignal::getRF()
{
    double w = 6.283185307179586/res.P;
    double c = 2.0*res.zeta*w;

    double *f = new double[n];

    for (int i = 0; i < n; i++)
    {
        f[i] = -ra[i] - c*rv[i];
    }

    return f;
}

double **EQSignal::getTHData()
{
    double **data = new double*[4];
    for (int i = 0; i < 4; i++) data[i] = new double[n];

    for (int i = 0; i < n; i++)
    {
        data[0][i] =   t[i];
        data[1][i] = acc[i];
        data[2][i] = vel[i];
        data[3][i] = disp[i];
    }
    return data;
}

void EQSignal::savetxt(const char *filename)
{
    ofstream out(filename, ios::out);

    out << "time\tacc\tvel\tdisp" << endl;

    for (int i = 0; i < n; i++)
    {
        out << t[i] << "\t" << acc[i] << "\t" << vel[i] << "\t" << disp[i] << endl;
    }

    out.close();

}

void EQSignal::savetxt(QString filename)
{
    this->savetxt(filename.toUtf8().data());
}

void EQSignal::savecsv(const char *filename)
{
    ofstream out(filename, ios::out);

    out << "time, acc, vel, disp" << endl;

    for (int i = 0; i < n; i++)
    {
        out << t[i] << ", " << acc[i] << ", " << vel[i] << ", " << disp[i] << endl;
    }

    out.close();

}

void EQSignal::savecsv(QString filename)
{
    this->savecsv(filename.toUtf8().data());
}

void EQSignal::savetxtsp(const char *filename)
{
    ofstream out(filename, ios::out);
    double ***data = this->getSPData();
    int np = sp->getNP();
    
    for (int i = 0; i < nsp; ++i)
    {
        out << "Period(zeta=" << zeta[i] <<")";
        out << "\tSPA\tSPV\tSPD\tSPE\tSPT";
        if (i<(nsp-1))
        {
            out << "\t";
        }
    }
    out << endl;

    for (int i = 0; i < np; ++i)
    {
        for (int j = 0; j < nsp; ++j)
        {
            for (int k = 0; k < 6; ++k)
            {
                out << data[j][k][i] <<"\t";
            }
        }
        out << endl;
    }

    out.close();

}

void EQSignal::savetxtsp(QString filename)
{
    this->savetxtsp(filename.toUtf8().data());
}

void EQSignal::savecsvsp(const char *filename)
{
    ofstream out(filename, ios::out);
    double ***data = this->getSPData();
    int np = sp->getNP();

    for (int i = 0; i < nsp; ++i)
    {
        out << "Period(zeta=" << zeta[i] <<")";
        out << ", SPA, SPV, SPD, SPE, SPT";
        if (i<(nsp-1))
        {
            out << ", ";
        }
    }

    out << endl;

    for (int i = 0; i < np; ++i)
    {
        for (int j = 0; j < nsp; ++j)
        {
            for (int k = 0; k < 6; ++k)
            {
                out << data[j][k][i] << ", ";
            }
        }
        out << endl;
    }

    out.close();
}

void EQSignal::savecsvsp(QString filename)
{
    this->savecsvsp(filename.toUtf8().data());
}

void EQSignal::savetxtsp(const char *filename, int indsp)
{
    ofstream out(filename, ios::out);

    Spectra *spi = sp+indsp;
    double **data = spi->getSPData();
    int np = spi->getNP();

    out << "Period\tSPA\tSPV\tSPD\tSPE\tSPT" << endl;

    for (int i = 0; i < np; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            out << data[i][j];
        }
        out << endl;
    }

    out.close();
}

void EQSignal::savetxtsp(QString filename, int indsp)
{
    this->savetxtsp(filename.toUtf8().data(),indsp);
}

void EQSignal::savecsvsp(const char *filename, int indsp)
{
    ofstream out(filename, ios::out);

    Spectra *spi = sp+indsp;
    double **data = spi->getSPData();
    int np = spi->getNP();

    out << "Period, SPA, SPV, SPD, SPE, SPT" << endl;

    for (int i = 0; i < np; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            out << data[i][j] << ", ";
        }
        out << endl;
    }

    out.close();
}

void EQSignal::savecsvsp(QString filename, int indsp)
{
    this->savecsvsp(filename.toUtf8().data(),indsp);
}

void EQSignal::setAcc(double *a)
{
	for (int i = 0; i < n; i++)
	{
		acc[i] = a[i];
	}
    confirm();
    a2vd();
}

void EQSignal::modifyAcc(double *a, int ind1, int ind2)
{
    int N = ind2<(n-1) ? ind2 : (n-1);
    for (int i = ind1; i <= N; i++)
    {
        acc[i] = a[i-ind1];
    }
    confirm();
    a2vd();
}


double **EQSignal::getSPData(int i)
{
    Spectra *spi = sp+i;
    return spi->getSPData();
}

void EQSignal::response(double zeta, double P, int method)
{
	switch (method)
	{
	case(0) :
		rfreq(acc, &n, &dt, &zeta, &P, ra, rv, rd);
		break;
	case(1) :
		rnmk(acc, &n, &dt, &zeta, &P, ra, rv, rd);
		break;
	case(2) :
		rmixed(acc, &n, &dt, &zeta, &P, ra, rv, rd);
		break;
	default:
		break;
	}

	fill(res.ku, n, 4.0*PI2 / (P*P));
    res.zeta = zeta;
    res.P = P;
}

void EQSignal::responseNL(double zeta, double P, int method, double *cp)
{
    int SM = method;
	rnl(acc, &n, &dt, &zeta, &P, ra, rv, rd, res.ku, &SM, cp);
    res.zeta = zeta;
    res.P = P;
    res.fy = cp[5];
    res.uy = cp[6];
}

double ***EQSignal::getSPData()
{
    double ***data = new double**[nsp];

    for (int i = 0; i < nsp; ++i)
    {
        data[i] = (sp+i)->getSPData();
    }

    return data;
}

double **EQSignal::getEnergy()
{
    double **e = new double*[5];
    for (int i = 0; i < 5; ++i) e[i] = new double[n];

    double *Ek = e[0];
    double *Ez = e[1];
    double *Es = e[2];
    double *Eh = e[3];
    double *Ein = e[4];

	energy(acc, &n, &dt, &(res.zeta), &(res.P), ra, rv, rd, Ek, Ez, Es, Eh, Ein, res.ku);

    return e;
}

void EQSignal::reallocateTH()
{
    
    if (0 != t)    delete [] t;
    if (0 != acc0) delete [] acc0;
    if (0 != acc)  delete [] acc;
    if (0 != vel)  delete [] vel;
    if (0 != disp) delete [] disp;
    if (0 != Ia)   delete [] Ia;
    if (0 != Iv)   delete [] Iv;
    if (0 != Id)   delete [] Id;
    if (0 != ta)   delete [] ta;
    if (0 != tv)   delete [] tv;
    if (0 != td)   delete [] td;
    if (0 != ra)   delete [] ra;
    if (0 != rv)   delete [] rv;
    if (0 != rd)   delete [] rd;
    if (0 != af)   delete [] af;
    if (0 != freqs)   delete [] freqs;
    if (0 != ampf)   delete [] ampf;
    if (0 != angf)   delete [] angf;
    if (0 != dangf)   delete [] dangf;
    if (0 != fpsd)   delete [] fpsd;
    if (0 != psd)   delete [] psd;

    t    = linspace(0.0,dt*(n-1),n);
    acc0 = new double[n];
    acc  = new double[n];
    vel  = new double[n];
    disp = new double[n];
    Ia = new double[n];
    Iv = new double[n];
    Id = new double[n];
    ta = new double[n];
    tv = new double[n];
    td = new double[n];
    ra = new double[n];
    rv = new double[n];
    rd = new double[n];

	res.ra = ra;
	res.rv = rv;
	res.rd = rd;
	res.ku = new double[n];

    nfft = nextpow2(n);
    npsd = nfft/NS;
    af = new std::complex<double>[nfft];
    double fn = 0.5/dt;
    freqs = linspace(0.0,fn,nfft/2+1);
    ampf = zeros(nfft/2+1);
    angf = zeros(nfft/2+1);
    dangf = zeros(nfft/2+1);
    fpsd = linspace(0.0,fn,npsd/2+1);
    psd = zeros(npsd/2+1);
}
